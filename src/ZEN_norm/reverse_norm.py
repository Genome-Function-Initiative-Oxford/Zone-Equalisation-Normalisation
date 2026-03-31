import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from .chrom_analysis import ChromAnalysisExtended

######################################
# Main class for reverse normalisation
######################################

class ReverseNorm(ChromAnalysisExtended):
    __slots__ = ("divisors", "contain_floats")
    def __init__(self, bigwig_paths = [], sample_names = [], blacklist = "", n_cores = 1, analysis_name = "Analysis", verbose = 1):
        """
        Reverse linear normalisation method applied to bigWigs.

        params:
            bigwig_paths:        List of paths of bigWig and/or wig files of interest. This takes
                                 priority over bigwig_directory.
            sample_names:        Optionally set as a list of custom names corresponding to each file.
                                 e.g. 'cntrl_s1.bw' and 'flt3_inhibitor_s1.bw' could be set as 
                                 ["Control Sample", "Treated Sample"].
                                 This will be converted to a dictionary mapping original file names to 
                                 the provided custom names.
                                 e.g. accessing sample_names would return 
                                 {"cntrl_s1": "Control Sample", "flt3_inhibitor_s1": "Treated Sample"}.
            blacklist:           File path to blacklist file with chromosome coordinates to exclude.
            n_cores:             The number of cores / CPUs to use if using multiprocessing.
            analysis_name:       Custom name of folder to save results to. By default this will be set 
                                 to "Analysis".
            verbose:             Set as an integer greater than 0 to display progress messages.
        """

        # Initialise the parent class
        super().__init__(bigwig_paths = bigwig_paths,
                         allow_bams = False, # Only allow bigWigs or wig files to be input
                         sample_names = sample_names,
                         blacklist = blacklist,
                         n_cores = n_cores,
                         analysis_name = analysis_name,
                         verbose = verbose)

        # Add output folders
        self.output_directories["bigwig"] = os.path.join(self.output_directories["results"], 
                                                         "BigWigs", "Reverse_Normalised")

        # Use to hold divisors per sample
        self.divisors = {}
        # Record whether each sample contains non-integer signal
        self.contain_floats = {}

    def signalReader(self, bw_idx, chromosome, start_idx = 0, end_idx = -1, pad_end = False,
                     dtype = np.float32, verbose = 0):
        """
        Read bigWig signal into a numpy array.

        params:
            bw_idx:      Index of sample to read signal for.
            chromosome:  Chromosome to read signal for.
            start_idx:   Start base-pair position (zero indexed).
            end_idx:     End base-pair position (zero indexed).
            pad_end:     Whether to add zeros to end to ensure array size is consistent across 
                         samples for a chromosome.
            dtype:       Data type to store the signal as.
            verbose:     Set as a number > 0 to print progress.
        """

        sample_name = list(self.sample_names.keys())[bw_idx]

        # Open the bigWig for the sample
        bigwig_file = self.bigwig_paths[bw_idx]
        
        return self.extractChunkSignal(bw_idx = bw_idx,
                                       bigwig_file = bigwig_file,
                                       sample_name = sample_name,
                                       chromosome = chromosome,
                                       start_idx = start_idx,
                                       end_idx = end_idx,
                                       pad_end = pad_end,
                                       dtype = dtype,
                                       verbose = verbose)
    
    def estimateFragmentSize(self, bw_idx, chromosome, start_idx = 0, end_idx = -1):
        """
        Get the smallest non-zero value from signal as a predictor of the signal magnitude 
        created by a signal read.

        params:
            bw_idx:     Sample ID to estimate fragment size for.
            chromosome: Chromosome to estimate fragment size for.
            start_idx:  Start base-pair position (zero indexed).
            end_idx:    End base-pair position (zero indexed).
        """

        # Assume all signal contains integers
        found_float = False

        # Get signal within range
        signal = self.signalReader(bw_idx = bw_idx, chromosome = chromosome, 
                                   start_idx = start_idx, end_idx = end_idx)
        # Drop zeros
        signal = signal[signal != 0]

        if len(signal) > 0:
            # Check the signal's sign
            if self.mean_genome_signals[bw_idx] < 0:
                # Signal is negative, so get the largest non-zero value
                fragment_estimate = max(signal)
            else:
                # Signal is positive, so get the smallest non-zero value
                fragment_estimate = min(signal)

            # Check if any of the signal is non-integer
            if not fragment_estimate.is_integer():
                # Fragment estimate was non-integer
                found_float = True
            else:
                # Check rest of signal
                signal = np.unique(signal)
                if not np.all(signal == np.floor(signal)):
                    found_float = True

        else:
            # All signal was zero, so cannot determine fragment size
            fragment_estimate = -np.inf

        return bw_idx, fragment_estimate, found_float

    def estimateFragmentSizes(self, chunk_size, chromosomes):
        """
        Get the smallest non-zero values from sample signals as predictors of the signal 
        magnitude created by each signal read.

        params:
            chunk_size:  Set as a non-negative number to estimate the fragment size by breaking 
                         the signal into chunks of this size, or as -1 to use whole chromosomes 
                         in one go.
            chromosomes: List of chromosomes to estimate fragment sizes for.
        """

        if self.verbose > 0:
            n_chroms = len(chromosomes)
            n_samples = len(self.sample_ids)
            print(f"Estimating fragment sizes from {n_chroms} chromosome{'s' if n_chroms != 1 else ''} "
                  f"for {n_samples} sample{'s' if n_samples != 1 else ''}")

        # Check whether to read chromosomes into smaller chunks
        use_whole_chroms = (chunk_size == -1)

        if len(chromosomes) == 0:
            chromosomes = self.sample_chromosomes
        else:
            missing_chroms = []

            for chrom in chromosomes:
                if chrom not in self.sample_chromosomes:
                    missing_chroms.append(chrom)

            n_missing = len(missing_chroms)

            if n_missing > 0:
                raise ValueError(f"{n_missing} chromosome{'s' if n_chroms != 1 else ''} were not found: "
                                 f'"{", ".join(missing_chroms)}"')

        # Record absolute values of smallest genome-wide fragment estimate per sample
        abs_frag_estimates = dict(zip(self.sample_ids, [np.inf] * len(self.sample_ids)))
        # Record if sample signal contains any non-integer values
        contain_floats = dict(zip(self.sample_ids, [False] * len(self.sample_ids)))

        if self.n_cores > 1:
            with ProcessPoolExecutor(self.n_cores) as executor:
                processes = []

                for chrom in chromosomes:
                    chrom_size = self.chrom_sizes[chrom]

                    if use_whole_chroms:
                        # Read the whole chromosome at once
                        chunk_size = chrom_size

                    for chunk_start in range(0, chrom_size, chunk_size):
                        chunk_end = min(chrom_size, chunk_start + chunk_size)

                        for bw_idx in self.sample_chromosomes[chrom]:
                            processes.append(executor.submit(self.estimateFragmentSize,
                                                             bw_idx = bw_idx,
                                                             chromosome = chrom,
                                                             start_idx = chunk_start,
                                                             end_idx = chunk_end))
                        chunk_start = chunk_end
                        chunk_end = min(chrom_size, chunk_start + chunk_size)
                        
                for process in as_completed(processes):
                    # Update the sample's fragment estimate as the smallest absolute value
                    bw_idx, frag_estimate, found_float = process.result()
                    abs_frag_estimates[bw_idx] = min(abs_frag_estimates[bw_idx], abs(frag_estimate))

                    if found_float:
                        contain_floats[bw_idx] = found_float

            if self.checkParallelErrors(processes):
                raise RuntimeError("estimateFragmentSizes failed to complete. To debug, see trace above.")
            
        else:
            for chrom in chromosomes:
                chrom_size = self.chrom_sizes[chrom]

                if use_whole_chroms:
                    chunk_size = chrom_size

                for chunk_start in range(0, chrom_size, chunk_size):
                    chunk_end = min(chrom_size, chunk_start + chunk_size)

                    for bw_idx in self.sample_chromosomes[chrom]:
                        bw_idx, frag_estimate, found_float = self.estimateFragmentSize(bw_idx = bw_idx,
                                                                                       chromosome = chrom,
                                                                                       start_idx = chunk_start,
                                                                                       end_idx = chunk_end)
                        abs_frag_estimates[bw_idx] = min(abs_frag_estimates[bw_idx], abs(frag_estimate))

                        if found_float:
                            contain_floats[bw_idx] = found_float

                    chunk_start = chunk_end
                    chunk_end = min(chrom_size, chunk_start + chunk_size)

        return abs_frag_estimates, contain_floats

    def calculateDivisors(self, chunk_size = 10000000, chromosomes = []):
        """
        Set divisors as the absolute fragment size estimates.

        params:
            chunk_size:  Set as a non-negative number to estimate the fragment size by breaking 
                         the signal into chunks of this size, or as -1 to use whole chromosomes 
                         in one go.
            chromosomes: List of chromosome(s) to predict fragment size from.
        """

        calculate_frags = True

        if (len(self.divisors) > 0) & (len(self.contain_floats) > 0):
            # Check for all values to divide by per samples
            if len(np.setdiff1d(list(self.divisors.keys()), self.sample_ids)) == 0:
                # Check if each sample has been classified as containing integers or floats
                if len(np.setdiff1d(list(self.contain_floats.keys()), self.sample_ids)) == 0:
                    # Found expected values
                    calculate_frags = False

        if calculate_frags:
            # Clear any previously set values
            self.divisors = {}
            self.contain_floats = {}
            # Calculate the absolute fragment size estimates to use as divisors
            self.divisors, self.contain_floats = self.estimateFragmentSizes(chunk_size, chromosomes)

    def testNorm(self, chunk_size = 10000000, chromosomes = []):
        """
        Test and classify samples based on whether each was already normalised.
        
        params:
            chunk_size:  Set as a non-negative number to estimate the fragment size by breaking 
                         the signal into chunks of this size, or as -1 to use whole chromosomes 
                         in one go.
            chromosomes: List of chromosome(s) to predict fragment size from.

        returns:
            norm_status: Dictinary classifying samples into unnormalised, confirmed non-linear 
                         normalised and linear / unconfirmed non-linear normalised
        """

        custom_sample_names = np.array(list(self.sample_names.values()))

        # Set divisors as the absolute fragment size estimates
        self.calculateDivisors(chunk_size, chromosomes)

        # Classify samples as either unnormalised or normalised
        raw_samples = []
        non_linear_samples = []
        norm_samples = []

        for bw_idx in self.sample_ids:
            # Get the sample specific divisor
            div = self.divisors[bw_idx]
            is_float = self.contain_floats[bw_idx]
            sample_name = custom_sample_names[bw_idx]

            if (div == 1) & (not is_float):
                raw_samples.append(sample_name)

            elif div == 1:
                non_linear_samples.append(sample_name)

            else:
                norm_samples.append(sample_name)

        # Classify samples into unnormalised, confirmed non-linear normalised and 
        # linear / unconfirmed non-linear normalised
        n_samples = len(custom_sample_names)
        n_raw = len(raw_samples)
        n_non_linear = len(non_linear_samples)
        n_norm = len(norm_samples)

        norm_status = {}

        if n_raw > 0:
            if self.verbose > 0:
                if n_raw == n_samples:
                    print("No samples were pre-normalised")
                else:
                    message = f"{n_raw} sample{'s' if n_raw != 1 else ''} are pre-normalised: "
                    message += ", ".join(raw_samples)
                    print(message)

            norm_status["Unnormalised"] = raw_samples

        if n_non_linear > 0:
            if self.verbose > 0:
                if n_non_linear == n_samples:
                    print("All samples were transformed by a non-linear normalisation method")
                else:
                    message = f"{n_non_linear} sample{'s' if n_non_linear != 1 else ''} were "
                    message += "transformed by a non-linear normalisation method: "
                    message += ", ".join(non_linear_samples)
                    print(message)

            norm_status["Non_Linear"] = non_linear_samples

        if n_norm > 0:
            if self.verbose > 0:
                if n_norm == n_samples:
                    print("All samples were pre-normalised")
                else:
                    message = f"{n_norm} sample{'s' if n_norm != 1 else ''} were "
                    message += "transformed by a normalisation method: "
                    message += ", ".join(norm_samples)
                    print(message)

            norm_status["Normalised"] = norm_samples

        return norm_status

    # Updated to reverse normalise the signal
    def openDefaultSignal(self, bw_idx, chromosome, signal_type = ""):
        """
        Open the raw signal and reverse the normalisation through division.

        params:
            bw_idx:      Sample ID to open signal for.
            chromosome:  Chromosome to open signal for.
            signal_type: Parameter is redundant for this class.
        """

        # Get the value to divide sample signal by
        div = self.divisors[bw_idx]

        # Open original signal
        signal = self.signalReader(bw_idx = bw_idx, chromosome = chromosome, start_idx = 0, end_idx = -1)
        # Reverse linear scaling
        signal = signal / div

        # Check floating point errors are within one decimal place
        signal_unique = np.sort(np.unique(signal[signal != 0,]))

        if np.any(np.abs(signal_unique - np.round(signal_unique)) > 0.1):
            sample_name = list(self.sample_names.values())[bw_idx]
            print(f"Warning: Floating point error was found within one decimal place for "
                  f"{sample_name} {chromosome}")

        signal = np.round(signal)

        return signal

    def runReversal(self, chunk_size = 10000000, chromosomes = [], replace_existing = False):
        """
        Run reverse normalisation of all sample bigWigs.
        
        params:
            chunk_size:       Set as a non-negative number to estimate the fragment size by breaking 
                              the signal into chunks of this size, or as -1 to use whole chromosomes 
                              in one go.
            chromosomes:      List of chromosome(s) to predict fragment size from. This differs
                              from the chromosomes that will be included in the reverse 
                              normalised bigWigs.
            replace_existing: Whether to overwrite previously created files.
        """

        custom_sample_names = np.array(list(self.sample_names.values()))

        if replace_existing:
            incomplete_bw_idxs = self.sample_ids

        else:
            # Find indexes of samples that have not yet been reverse normalised
            incomplete_bw_idxs = self.findCurrentFiles(directories = [self.output_directories["bigwig"]],
                                                       file_names_regex = "_reverse_norm.bw",
                                                       replace_existing = replace_existing)

            if len(incomplete_bw_idxs) == 0:
                if self.verbose > 0:
                    print("Reverse normalised bigWigs already created for all samples")
                
                return None

        # Set divisors as the absolute fragment size estimates
        self.calculateDivisors(chunk_size, chromosomes)

        if self.n_cores > 1:
            processes = []
            executor = ProcessPoolExecutor(self.n_cores)

        reversed_ids = []

        for bw_idx in incomplete_bw_idxs:
            # Get the sample specific divisor
            div = self.divisors[bw_idx]
            is_float = self.contain_floats[bw_idx]
            sample_name = custom_sample_names[bw_idx]

            if (div == 1) & (not is_float):
                if self.verbose > 0:
                    print(f"{sample_name} is not normalised as it has a fragment estimate of 1 and all "
                          f"values are integers ")

            elif div == 1:
                if self.verbose > 0:
                    # Sample has non-integer values, yet fragment estimate is 1
                    print(f"Cannot reverse {sample_name} non-linear normalisation for {sample_name}")

            else:
                if self.verbose > 0:
                    print(f"{sample_name} was pre-normalised and has a divisor of {round(div, 3)}")

                reversed_ids.append(bw_idx)

                if self.n_cores > 1:
                    processes.append(executor.submit(self.saveBigWig,
                                                     bw_idx = bw_idx,
                                                     custom_sample_name = sample_name,
                                                     file_name = f"{sample_name}_reverse_norm.bw",
                                                     directory = self.output_directories["bigwig"]))
                else:
                    self.saveBigWig(bw_idx = bw_idx,
                                    custom_sample_name = sample_name,
                                    file_name = f"{sample_name}_reverse_norm.bw",
                                    directory = self.output_directories["bigwig"])
                    
        if len(reversed_ids) > 0:
            if self.n_cores > 1:
                if self.checkParallelErrors(processes):
                    raise RuntimeError("runReversal failed to complete. To debug, see trace above.")
                
            if self.verbose > 0:
                print("Finished saving bigWigs")
