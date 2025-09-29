import numpy as np
import pandas as pd
import pyBigWig
import os
import glob
import math
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import subprocess
from seaborn import color_palette
from concurrent.futures import ProcessPoolExecutor, as_completed

#############################################################################################
# Baseline classes for ZEN normalisation, reverse normalisation and normalisation comparision
#############################################################################################

class ChromAnalysisCore:
    # Define attributes in slots to reduce RAM usage
    __slots__ = ("sample_names", "chromosomes", "chrom_sizes", "output_directories", "n_cores", 
                 "analysis_name", "verbose")
    def __init__(self, n_cores = 1, analysis_name = "Analysis", verbose = 1):
        """
        Core parent class for analysing bigWig signal.

        params:
            n_cores:             The number of cores / CPUs to use if using multiprocessing.
            analysis_name:       Custom name of folder to save results to. By default this will
                                 be set to "Analysis".
            verbose:             Set as an integer greater than 0 to display progress messages.
        """

        # Set base parameters
        self.setNCores(n_cores)
        self.setAnalysisName(analysis_name)
        self.setVerbose(verbose)

        # Directories to save intermediate files and results
        temp_dir = os.path.join(self.analysis_name, "Temp")
        results_dir = os.path.join(self.analysis_name, "Results")

        self.output_directories = {"temp": temp_dir,
                                   "results": results_dir}

    def __str__(self):
        """ Format parameters in user readable way when print is called on an object """

        message = (f'{self.__class__.__name__} object for "{self.analysis_name}"\n'
                   f"   * Resources: {self.n_cores} cores\n")

        return message

    def setNCores(self, n_cores):
        try:
            self.n_cores = max(int(n_cores), 1)
        except:
            raise ValueError("n_cores must be a non-negative integer as this is the number of "
                             "CPUs to use for parallelisation.")

    def setAnalysisName(self, analysis_name):
        if (analysis_name is None) or (analysis_name == ""):
            # Set a name for the analysis results folder
            self.analysis_name = "Analysis"
        else:
            self.analysis_name = str(analysis_name)

    def setVerbose(self, verbose):
        try:
            self.verbose = int(verbose)
        except:
            self.verbose = 1

    @staticmethod
    def checkValidDirectory(directory):
        """
        Raises error if directory is missing or empty.
        """

        if not os.path.exists(directory):
            raise ValueError(f'Directory does not exist: "{directory}"')
        elif len(os.listdir(directory)) == 0:
            raise FileNotFoundError(f'Directory is empty: "{directory}"')

    @staticmethod
    def naturalSort(values, return_idxs = False):
        """
        Sort strings with numbers inside,
        e.g. ["text_1", "text_100", "text_2"] to ["text_1", "text_2", "text_100"].

        params:
            values:      List or array of strings.
            return_idxs: Set as True to also return the indexes.

        returns:
            sorted_values: Array of sorted strings.
        """

        values = np.array(values)
        sorting_keys = []

        for val in values:
            sorting_keys.append([int(x) if x.isdigit() else x for x in re.split(r'(\d+)', val)])

        sorting_idxs = sorted(range(len(sorting_keys)), key = lambda i: sorting_keys[i])
        sorted_values = values[sorting_idxs]

        if return_idxs:
            return sorted_values, sorting_idxs
        else:
            return sorted_values

    @staticmethod
    def checkParallelErrors(parallel):
        """
        Check for unreported errors caused by multithreading or multiprocessing.
        
        params:
            parallel: List of threads / processes to check for errors.
        """

        # Detect any exceptions raised by threads / processes
        exceptions = [p.exception() for p in parallel]
        raised_exceptions = [e for e in exceptions if e is not None]

        if not raised_exceptions:
            return False

        n_failed = len(raised_exceptions)
        print(f"\nError {n_failed} process{'es' if n_failed != 1 else ''} failed:")

        # Count number of times each exception type occured
        exception_counts = {}

        for e in raised_exceptions:
            key = f"{e.__class__.__name__}: {e}"
            if key in exception_counts:
                exception_counts[key] += 1
            else:
                exception_counts[key] = 1

        for e, count in exception_counts.items():
            print(f"{count} process{'es' if count != 1 else ''} raised: {e}")

        if any(e.__class__.__name__ == "BrokenProcessPool" for e in raised_exceptions):
            print("'BrokenProcessPool' can sometimes be caused by a lost connection or inadequate ",
                  "resources.",
                  "For example, requesting n_cores=10 when only 8 CPUs are avaliable to the system.\n")

        return True

    @staticmethod
    def openFiles(file_paths = np.array([]), directory = "", min_files = 1, verbose = 0):
        """
        Open all bigWig files in a directory.
        
        params:
            file_paths:     A list of paths to bigWig files of interest.
            directory:      The working directory containing the bigWig files of interest.
            min_files:      Throw an error if less files are found than this number.
            verbose:        Set as > 0 to print messages for debugging.
        """

        files = []

        file_extensions = [".bw", ".bigWig"]

        if len(directory) != 0:
            # Find the files within the directory
            file_paths = np.array([])

            for extension in file_extensions:
                file_paths = np.concatenate((file_paths, 
                                             glob.glob(os.path.join(directory, ("*" + extension)))))

            if len(file_paths) == 0:
                # Attempt to add working directory if cannot find files
                directory = os.path.join(os.getcwd(), directory)
                for extension in file_extensions:
                    file_paths = np.concatenate((file_paths,
                                                 glob.glob(os.path.join(directory, ("*" + extension)))))

            if len(file_paths) < max(min_files, 1):
                raise ValueError(f"{len(file_paths)} bigwig files found using extension"
                                 f"{'s ' if len(file_extensions) > 1 else ' '}"
                                 f"{', '.join(e for e in file_extensions)} "
                                 f'in directory "{directory}".\n'
                                 f"This is less than the minimum expected number of {min_files}.\n"
                                 f"You may alternatively set a list of files manually via file_paths.")
            
        elif len(file_paths) > 0:
            file_paths_filtered = np.array([])

            for extension in file_extensions:
                if isinstance(file_paths, str):
                    # If single file given as a string, check if it has the correct extension
                    if file_paths.endswith(extension):
                        file_paths_filtered = np.array([file_paths])
                
                else:
                    # If array of files given, keep only those with the correct extension
                    file_paths_filtered = np.concatenate((file_paths_filtered, 
                                                          file_paths[np.char.endswith(file_paths, 
                                                                                      extension)]))

            file_paths = file_paths_filtered

            if len(file_paths) < max(min_files, 1):
                raise ValueError(f"{len(file_paths)} bigWig files found in file_paths.\n"
                                 f"This is less than the minimum expected number of {min_files}.\n"
                                 f"Ensure this is set as a list of files with extension"
                                 f'{", ".join(e for e in file_extensions)}, '
                                 f"or alternatively set parameter directory as a directory containing "
                                 f"the files.")
            
        else:
            raise ValueError("Set parameter directory as a directory containing the files or "
                             "file_paths as a list of files to find.")

        if verbose > 0:
            print(f"Opening {len(file_paths)} bigWig files")

        for path in file_paths:
            if os.path.exists(path):
                try:
                    files.append(pyBigWig.open(path))
                except Exception as e:
                    print(e)
                    raise Exception(f'Could not open "{path}"')
            else:
                raise FileNotFoundError(f'Could not find "{path}"')

        return files, file_paths

    def orderChroms(self, chromosomes):
        """
        Sort chromosomes by numerical and non-numerical names.

        params:
            chromosomes: List of chromosomes to order.
        """

        # Search for a common prefix, e.g. for ["chr11", "chr7", "chrY", "chr2", "chrX"] the prefix 
        # is "chr"
        chrom_prefix = os.path.commonprefix(list(chromosomes))
        # Remove the prefix from each chromosome, e.g. ["11", "7", "Y", "2", "X"]
        chrom_postfixes = [chrom.replace(chrom_prefix, '', 1) for chrom in chromosomes]

        # Separate into numerical and non-numerical postfixes, e.g. [11, 7, 2] and ["Y", "X"]
        numeric_postfixes = []
        character_postfixes = []

        for postfix in chrom_postfixes:
            if postfix.isdigit():
                numeric_postfixes.append(int(postfix))
            else:
                character_postfixes.append(postfix)

        # Sort the postfixes, e.g. [2, 7, 11] and ["X", "Y"]
        numeric_postfixes = np.sort(numeric_postfixes)
        character_postfixes = np.sort(character_postfixes)

        # Create ordered array of chromosomes, e.g. ["chr2", "chr7", "chr11", "chrX", "chrY"]
        ordered_chromosomes = [chrom_prefix + str(num) for num in numeric_postfixes]
        ordered_chromosomes.extend([chrom_prefix + char for char in character_postfixes])
        ordered_chromosomes = np.array(ordered_chromosomes)

        return ordered_chromosomes

    def calculateChromSizes(self, bigwigs, chromosomes, get_sample_chroms = False):
        """
        Set a dictionary mapping chromosomes to number of base-pairs dervied from bigWigs. 
        Where chromosome size differs amoung samples, the largest size is used.

        params:
            bigwigs:           List of opened pyBigWig bigWigs.
            chromosomes:       List of chromosomes to get the sizes for.
            get_sample_chroms: Set as True to also create a dictionary of sample-specific 
                               chromosome sizes.
        """

        # Get the size of each chromosome
        chrom_sizes = {}

        if get_sample_chroms:
            sample_chromosomes = {}

        if self.verbose > 0:
            unmatched_chroms = []

        all_bw_chroms = [bw.chroms() for bw in bigwigs]

        for chrom in chromosomes:
            size_values = []

            # Iterate through bigwigs
            for bw_idx, bw_chroms in enumerate(all_bw_chroms):
                chrom_size = bw_chroms.get(chrom)

                if chrom_size is not None:
                    # Record the size of the chromosome saved within the bigWig
                    size_values.append(chrom_size)

                    if get_sample_chroms:
                        if bw_idx == 0:
                            # Record first index
                            sample_chromosomes[chrom] = [bw_idx]
                        else:
                            # Append further indexes
                            sample_chromosomes[chrom].append(bw_idx)
                
            size_values = np.unique(size_values)

            if len(size_values) > 1:
                size_values = np.sort(size_values)
                if self.verbose > 0:
                    unmatched_chroms.append(chrom)

            # Set the largest chromosome size
            chrom_sizes[str(chrom)] = int(size_values[-1])
                
        if (self.verbose > 0):
            if len(unmatched_chroms) > 0:
                print(f"Warning: {len(unmatched_chroms)} chromosome(s) were found to have differing "
                      f"sizes across samples: "
                      f'{", ".join(i for i in unmatched_chroms)}')

        if get_sample_chroms:
            return chrom_sizes, sample_chromosomes
        else:
            return chrom_sizes

    def extractChunkSignal(self, bigwig_file, chromosome, start_idx, end_idx, pad_end = False, 
                           dtype = np.float32, bw_idx = None, sample_name = "", 
                           return_file = False, verbose = 0):
        """
        Read signal between given coordinates from bigWig.

        params:
            bigwig_file: File path to a bigWig to get signal from.
            chromosome:  Chromosome to read signal from.
            start_idx:   Start base-pair position (zero indexed).
            end_idx:     End base-pair position (zero indexed).
            pad_end:     Whether to add additonal zeros to end if reading full chromosome and
                         sample's chromosome size is smaller than for other samples.
            dtype:       Data type to store the signal as.
            bw_idx:      Sample ID optionally given for printing.
            sample_name: Sample name optionally given for printing.
            return file: Set as True to return the bigwig file alongside the signal.
            verbose:     Set as a number > 0 to print progress.
        """

        if not os.path.exists(bigwig_file):
            raise ValueError(f'File "{bigwig_file}" does not exist')

        if sample_name:
            message_postfix = f" for {sample_name}"
        else:
            message_postfix = ""

        bigwig = pyBigWig.open(bigwig_file)

        if chromosome not in bigwig.chroms():
            raise ValueError(f"Chromosome {chromosome} was not found{message_postfix}")
        
        if (end_idx > 0) and (start_idx >= end_idx):
            if verbose > 0:
                print(f"End exceeded start{message_postfix} {chromosome}:{start_idx}-{end_idx}")
            return np.array([], dtype = dtype)
        
        elif start_idx >= self.chrom_sizes[chromosome]:
            if verbose > 0:
                print(f"Start exceeded chromosomes size of "
                      f"{self.chrom_sizes[chromosome]}{message_postfix} "
                      f"{chromosome}:{start_idx}-{end_idx}")
            return np.array([], dtype = dtype)

        bw_chrom_size = bigwig.chroms()[chromosome]
        array_end = end_idx

        if (end_idx < 0) or (end_idx > bw_chrom_size):
            # Cap size to fit chromosome
            end_idx = bw_chrom_size

            if pad_end and (chromosome in self.chrom_sizes):
                # Ensure array size is consistent across samples
                array_end = self.chrom_sizes[chromosome]
            else:
                array_end = end_idx

        if start_idx >= array_end:
            if verbose > 0:
                print(f"End exceeded start{message_postfix} {chromosome}:{start_idx}-{end_idx}")
            return np.array([], dtype = dtype)

        if verbose > 0:
            message = f"Reading {chromosome}:{start_idx}-{end_idx}{message_postfix}"

            if bw_idx is not None:
                message += f" ({bw_idx + 1}/{len(self.sample_names)})"
            
            print(message)

        # Create placeholder array of zeros
        chunk_signal = np.zeros(array_end - start_idx, dtype = dtype)
        # Extract signal for the chromosome chunk
        chunk_signal[:(end_idx - start_idx)] = np.array(bigwig.values(chromosome, start_idx, end_idx), 
                                                        dtype = dtype)
        # Ensure zeros are +0 not -0
        chunk_signal[(chunk_signal == 0) & np.signbit(chunk_signal)] = 0
        # Replace any nan values with zero
        chunk_signal[np.isnan(chunk_signal)] = 0

        if return_file:
            # Return file as well as signal as identifier for multiprocessing
            return bigwig_file, chunk_signal
        else:
            return chunk_signal
        
class ChromAnalysisExtended(ChromAnalysisCore):
    # Define attributes in slots to reduce RAM usage
    __slots__ = ("bam_paths", "bam_directory", "bigwig_paths", "bigwig_directory", "wig_paths",
                 "file_sample_names", "sample_chromosomes", "mean_genome_signals", "sample_ids", 
                 "blacklist")
    def __init__(self, bam_paths = [], bam_directory = "", bigwig_paths = [], bigwig_directory = "",
                 allow_bams = False, sample_names = [], chromosomes = [], blacklist = "", n_cores = 1, 
                 analysis_name = "Analysis", verbose = 1):
        """
        Build upon the parent class to enable signal to be read from either bigWigs, wigs or BAMs.

        params:
            bam_paths:           List of paths of BAM files of interest. This takes priority over 
                                 bam_directory, bigwig_paths and bigwig_directory.
            bam_directory:       Directory containing the BAM files of interest. This takes priority
                                 over bigwig_paths and bigwig_directory.
            bigwig_paths:        List of paths of bigWig and/or wig files of interest. This takes 
                                 priority over bigwig_directory.
            bigwig_directory:    Directory containing the bigWig and/or wig files of interest. This
                                 is required if neither bam_paths, bam_directory or bigwig_paths are
                                 set.
            allow_bams:          Setting this as True allows BAMs to be supported as an input. 
                                 Designed for setting up subclasses.
            sample_names:        Optionally set as a list of custom names corresponding to each file.
                                 e.g. 'cntrl_s1.bw' and 'flt3_inhibitor_s1.bw' could be set as 
                                 ["Control Sample", "Treated Sample"].
                                 This will be converted to a dictionary mapping original file names
                                 to the provided custom names.
                                 e.g. accessing sample_names would return 
                                 {"cntrl_s1": "Control Sample", "flt3_inhibitor_s1": "Treated Sample"}.
            chromosomes:         List of chromosomes to run analysis on.
            blacklist:           File path to blacklist file with chromosome coordinates to exclude.
            n_cores:             The number of cores / CPUs to use if using multiprocessing.
            analysis_name:       Custom name of folder to save results to. By default this will be set 
                                 to "Analysis".
            verbose:             Set as an integer greater than 0 to display progress messages.
        """

        # Initialise the parent class
        super().__init__(n_cores = n_cores,
                         analysis_name = analysis_name,
                         verbose = verbose)

        # Set the directory or paths of input files and the sample names
        self.setSamples(bam_paths, str(bam_directory), bigwig_paths, str(bigwig_directory), 
                        allow_bams, sample_names)
        self.setBlacklist(blacklist)

        # Can only run straight away if bigWigs were input
        if not allow_bams:
            self.setChromAttributes(chromosomes)

    def __str__(self):
        """ 
        Format parameters in user readable way when print is called on an object.
        """

        n_chroms = len(self.chromosomes)
        if n_chroms == 0:
            chrom_msg = "Chromosomes: Not yet set"
        elif n_chroms == 1:
            chrom_msg = f"1 chromosome: {self.chromosomes[0]}"
        else:
            chrom_msg = f"{n_chroms} chromosome(s): {', '.join(self.chromosomes)}"

        message = (f'{self.__class__.__name__} object for "{self.analysis_name}"\n'
                   f"   * {chrom_msg}\n"
                   f"   * Number of samples: {len(self.sample_ids)}\n"
                   f"   * Sample names: {', '.join(self.sample_names.values())}\n"
                   f"   * Resources: {self.n_cores} cores\n")

        return message

    def setSamples(self, bam_paths, bam_directory, bigwig_paths, bigwig_directory, allow_bams, 
                   sample_names):
        """
        Check if BAM files (if allowed) or wigs/bigWigs were provided as either a path or directory.
        Then keep either a directory or path (path > directory) and extract the sample names for 
        the files.
        """

        if allow_bams:
            self.bam_paths = np.array([], dtype = str)
            self.bam_directory = ""

        self.bigwig_paths = np.array([], dtype = str)
        self.bigwig_directory = ""
        self.wig_paths = np.array([], dtype = str)

        # Paths take priority over a directory
        if (allow_bams and (len(bam_paths) > 0)):
            if isinstance(bam_paths, str) and os.path.isdir(bam_paths):
                # Directory was given, so use this rather than list of file paths
                bam_directory = bam_paths
            else:
                self.bam_paths = np.array(list(bam_paths))

        elif (len(bigwig_paths) > 0):
            if isinstance(bigwig_paths, str) and os.path.isdir(bigwig_paths):
                bigwig_directory = bigwig_paths
            else:
                self.bigwig_paths = np.array(list(bigwig_paths))

        # BAM directory takes priority over bigWig directory
        if allow_bams and bam_directory and (len(self.bam_paths) == 0):
            # Check that directory exists and is not empty
            self.checkValidDirectory(bam_directory)

            # Get file paths of BAM files in the inputted directory
            self.bam_paths = np.array(glob.glob(os.path.join(bam_directory, "*.bam")))
            self.bam_directory = bam_directory

            if len(self.bam_paths) == 0:
                raise ValueError(f'No BAM files found in directory "{bam_directory}"')

        elif bigwig_directory and (len(self.bigwig_paths) == 0):
            self.checkValidDirectory(bigwig_directory)
            self.bigwig_directory = bigwig_directory
            # Check if any wig files were given
            wig_files = glob.glob(os.path.join(self.bigwig_directory, "*.wig"))

            if len(self.wig_paths) > 0:
                # Set wig files to convert to bigWig
                self.wig_paths = wig_files

            self.bigwig_paths = np.concatenate((glob.glob(os.path.join(self.bigwig_directory, "*.bw")),
                                                glob.glob(os.path.join(self.bigwig_directory, "*.bigWig"))))

            if (len(self.bigwig_paths) == 0) and (len(self.wig_paths) == 0):
                raise ValueError(f'No bigWig or wig files found in directory "{self.bigwig_directory}"')

        if (len(self.bigwig_paths) == 0) and (len(self.wig_paths) == 0):
            if allow_bams:
                if (len(self.bam_paths) == 0):
                    raise ValueError("Either BAM, bigWig or wig files need to be given by setting one "
                                     "of: bam_paths, bam_directory, bigwig_paths, bigwig_directory")
            else:
                raise ValueError("Either bigWig or wig files need to be given by setting bigwig_paths "
                                 "or bigwig_directory")

        if allow_bams and (len(self.bam_paths) > 0):
            # Derive sample names from BAM files
            full_sample_names = []
            missing_bam_indexes = []

            for path in self.bam_paths:
                # Extract sample name from file path
                path = path.split(os.sep)
                sample_name = path[-1].replace(".bam", "")
                full_sample_names.append(sample_name)
                # Check for BAM index
                bam_index_path = os.path.join(os.sep.join(path[:-1]), f"{sample_name}.bam.bai")

                if not os.path.exists(bam_index_path):
                    missing_bam_indexes.append(sample_name)
            
            if missing_bam_indexes:
                raise ValueError(f"BAM indexes missing for: {', '.join(missing_bam_indexes)}")

        else:
            # Find all values in file path with a wig or bigWig extension
            wig_files = self.bigwig_paths[np.char.endswith(self.bigwig_paths, ".wig")]
            bigwig_files = np.concatenate((self.bigwig_paths[np.char.endswith(self.bigwig_paths, ".bw")],
                                           self.bigwig_paths[np.char.endswith(self.bigwig_paths, ".bigWig")]))
            
            if (len(wig_files) + len(bigwig_files)) != len(self.bigwig_paths):
                raise ValueError(f"bigwig_paths contains files that are neither wig nor bigWig: "
                                 f'"{", ".join(self.bigwig_paths)}"')

            # Derive sample names from wig and bigWig files
            full_sample_names = []

            if len(wig_files) > 0:
                # Set wig files to convert to bigWig
                self.wig_paths = wig_files

                full_sample_names = [path.split(os.sep)[-1].replace(".wig", "") for 
                                     path in wig_files]

            if len(bigwig_files) > 0:
                full_sample_names.extend([path.split(os.sep)[-1].replace(".bw", "").replace(".bigWig", "") for 
                                          path in bigwig_files])

        # Sort names alphabetically and numerically
        full_sample_names, sort_idxs = self.naturalSort(full_sample_names, return_idxs = True)
        full_sample_names = [str(s) for s in full_sample_names]
        # Replace any underscores or dots as these are excluded from file names
        self.file_sample_names = [s.replace("_", "-").replace(".", "-") for s in full_sample_names]

        if len(full_sample_names) != len(np.unique(self.file_sample_names)):
            raise ValueError("Sample are named too similar as they need to be distinguishable without "
                             "_ and . characters. Please rename the sample files.")

        # Sort file paths to match order
        if len(self.bigwig_paths) > 0:
            self.bigwig_paths = self.bigwig_paths[sort_idxs]
        if allow_bams:
            if len(self.bam_paths) > 0:
                self.bam_paths = self.bam_paths[sort_idxs]

        if len(sample_names) == len(full_sample_names):
            # Match custom names with file names
            self.sample_names = dict(zip(full_sample_names, sample_names[sort_idxs]))
            if self.verbose > 1:
                print(f"Custom sample names were mapped as follows: {self.sample_names}")

        elif len(sample_names) > 0:
            raise ValueError(f"Could not match custom sample names with files found.\n"
                             f"{len(full_sample_names)} files were found, but sample_names "
                             f"is of length {len(sample_names)}.")
        
        else:
            # Use the full file names
            self.sample_names = dict(zip(full_sample_names, full_sample_names))

        # Generate numerical ID / index for each sample
        # e.g. [0,1,2] for samples named ["Patient_1", "Patient_2", "Control"]
        self.sample_ids = np.arange(len(self.sample_names), dtype = np.uint16)

    def setBlacklist(self, blacklist):
        if (blacklist is None) or (blacklist == ""):
            self.blacklist = None
        else:
            # Check blacklist is valid file
            self.blacklist = str(blacklist)
            if not os.path.exists(blacklist):
                raise ValueError(f'Could not find blacklist file "{blacklist}"')

    def checkChromosomes(self, bigwigs, chromosomes):
        """ 
        Check that all chromosomes set by the user to analyse are valid and present within 
        the bigWig files.

        params:
            bigwigs: List of pyBigWig bigwig files.
            chromosomes: List of chromosomes.

        returns:
            chrom_lengths: Dictionary of valid chromosomes and their lengths.
        """

        if self.verbose > 0:
            if len(chromosomes) > 0:
                print("Checking provided chromosomes are valid")
            else:
                print("Setting chromosomes as autosomal and sex chromosomes")

        # Get the names of all chromosomes within the bigWigs
        all_found_chroms = set()

        for i in range(len(bigwigs)):
            all_found_chroms.update(set(bigwigs[i].chroms()))

        if len(chromosomes) > 0:
            chromosomes = set(chromosomes)
            # Check if any of the provided chromosomes to query are absent
            missing_chroms = np.array(list(chromosomes.difference(all_found_chroms)))
            # Record which were successfully found
            query_chroms = chromosomes & all_found_chroms
            raise_error = False

            if len(missing_chroms) > 0:
                if np.all(np.char.startswith(missing_chroms, "chr") == False):
                    # Add prefix to names and test to see if chromosome was named incorrectly
                    missing_chroms_renamed = set(["chr" + str(chrom) for chrom in missing_chroms])
                    
                    if len(np.array(list(missing_chroms_renamed.difference(all_found_chroms)))) == 0:
                        # Update the chromosomes to query if successfully found matching names for 
                        # them all
                        query_chroms.update(missing_chroms_renamed)
                    else:
                        raise_error = True
                else:
                    raise_error = True

            if raise_error:
                all_found_chroms = np.array(list(all_found_chroms))
                all_found_chroms.sort()

                raise ValueError(f"The following chromosomes were not found in the bigWig files: "
                                 f"{', '.join(missing_chroms)}.\n"
                                 f"Ensure they are named correctly and that you have provided the "
                                 f"right files.\n"
                                 f"The names of valid chromosomes found across the files are:"
                                 f"\n{', '.join(all_found_chroms)}")
        else:
            # Find all autosomal and sex chromosomes
            chromosomes = np.array([chrom for chrom in all_found_chroms if re.match(r'\Achr(\d+|o|w|x|y|z)$', 
                                                                                    chrom, re.IGNORECASE)])
            if len(chromosomes) == 0:
                # Try again without the chr prefix
                chromosomes = np.array([chrom for chrom in all_found_chroms if re.match(r'\A(\d+|o|w|x|y|z)$', 
                                                                                        chrom, re.IGNORECASE)])

            query_chroms = chromosomes
            
        # Reorder the chromosomes by name
        ordered_chromosomes = self.orderChroms(chromosomes = query_chroms)

        # Get size of each chromosome and determine what chromosomes each sample has
        chrom_sizes, sample_chromosomes = self.calculateChromSizes(bigwigs = bigwigs, 
                                                                   chromosomes = ordered_chromosomes,
                                                                   get_sample_chroms = True)

        self.chrom_sizes = chrom_sizes
        self.sample_chromosomes = sample_chromosomes

    def calculateMeanGenomeSignal(self, bigwigs):
        """
        Calculate the mean signal across the genome per bigWig.

        params:
            bigwigs: List of opened bigWig files to extract stats from.
        """

        # Global mean signal per sample
        mean_genome_signals = np.zeros(len(bigwigs))

        for bw_idx in range(len(bigwigs)):
            # Calculate the average signal across all chromosomes in the bigWig
            bigwig = bigwigs[bw_idx]
            sum_data = bigwig.header()["sumData"]
            n_bases_covered = bigwig.header()["nBasesCovered"]

            if (sum_data == 0) or (n_bases_covered == 0):
                mean_genome_signals[bw_idx] = 0
            else:
                mean_genome_signals[bw_idx] = bigwig.header()["sumData"] / bigwig.header()["nBasesCovered"]

        return mean_genome_signals

    def setChromAttributes(self, chromosomes):
        """
        Set chromosome sizes and mean genome signal per bigWig for quality filtering.

        params:
            chromosomes: List of chromosomes.
        """

        if len(self.bigwig_paths) > 0:
            # Open each bigWig file in the list of paths
            bigwigs, self.bigwig_paths = self.openFiles(file_paths = self.bigwig_paths, 
                                                        min_files = 1, verbose = self.verbose)
        elif len(self.bigwig_directory) > 0:
            # Open each bigWig file in the directory
            bigwigs, self.bigwig_paths = self.openFiles(directory = self.bigwig_directory, 
                                                        min_files = 1, verbose = self.verbose)
        else:
            raise ValueError("Could not determine which bigWig files to run analysis on. "
                             "Please specify bigWigs of interest by setting either:\n"
                             "(a) bigwig_paths as a list/array of bigWig file paths\n"
                             "(b) bigwig_directory as a path to a folder containing the bigWig files")
        
        # Convert chromosomes to a set to prevent duplicates
        if isinstance(chromosomes, str):
            self.chromosomes = {chromosomes}
        else:
            self.chromosomes = set(chromosomes)

        # Set chromosome sizes
        self.checkChromosomes(bigwigs, chromosomes)
        self.chromosomes = np.array(list(self.chrom_sizes.keys()))
        # Calculate the global mean signal per sample
        self.mean_genome_signals = self.calculateMeanGenomeSignal(bigwigs)

    def updateSampleName(self, full_name, custom_name):
        """ Set a custom name for a sample.
        
        params:
            full_name:   The sample name matching the bigWig,
                         e.g. "b_cell_rep_1_1.bw" becomes "b_cell_rep_1_1".
            custom_name: Alternative name for the sample.
        """
        if full_name in self.sample_names.keys():
            self.sample_names[full_name] = custom_name
        else:
            raise IndexError(f"Sample {full_name} not found in sample_names")

        if self.verbose > 1:
            print(f"{full_name} was given the custom name {custom_name}")


    def regexFindSamples(self, regex, second_regex = "", ignore_case = False, 
                         return_custom_names = True, exclude_samples = np.array([])):
        """
        Find all sample names that match regular expressions.

        params:
            regex:               Regex to match to sample names. For example, "b_cell" will match 
                                 "treated_b_cell_1" and "b_cell_control" but not "t_cell_control_1".
            second_regex:        Optional additional regex to match to sample names. Setting this 
                                 will enable pairs of sample names to be returned formed comparing 
                                 the first and second regex. If this is set as the same value as
                                 regex, all paired combination of samples that match regex are returned.
            ignore_case:         Set as True to make the search case-insensitive.
            return_custom_names: Set as True to return custom set sample names, or False to return 
                                 full sample names.
            exclude_samples:     Optionally set list of sample names to ignore.
        """

        if ignore_case:
            # Make case-insensitive
            re_case = re.IGNORECASE
        else:
            re_case = 0

        # Get sample name variants
        full_sample_names = np.array(list(self.sample_names.keys()))
        custom_sample_names = np.array(list(self.sample_names.values()))
        exclude_samples = np.array(exclude_samples)
        n_exclude = 0

        if len(exclude_samples) > 0:
            # Find indexes of the sample names to ignore
            exclude_idxs = np.union1d(np.where(np.isin(full_sample_names, exclude_samples))[0],
                                      np.where(np.isin(custom_sample_names, exclude_samples))[0])
            n_exclude = len(exclude_idxs)

            if n_exclude > 0:
                # Remove names of samples
                full_sample_names = np.delete(full_sample_names, exclude_idxs)
                custom_sample_names = np.delete(custom_sample_names, exclude_idxs)

        # Record sample names that match the first regex
        regex_samples = []
        unique_second_regex = False

        if second_regex != "":
            create_pairs = True
            if regex != second_regex:
                # If second unique regex is given, record sample names that match it
                second_regex_samples = []
                unique_second_regex = True
        else:
            create_pairs = False
        
        # Check each sample name to see if it falls within either regular expression
        for full_name, custom_name in zip(full_sample_names, custom_sample_names):
            check_1 = re.search(regex, full_name, re_case)
            check_2 = re.search(regex, custom_name, re_case)
            if check_1 or check_2:
                if return_custom_names:
                    regex_samples.append(custom_name)
                else:
                    regex_samples.append(full_name)

            elif unique_second_regex:
                check_1 = re.search(second_regex, full_name, re_case)
                check_2 = re.search(second_regex, custom_name, re_case)
                if check_1 or check_2:
                    if return_custom_names:
                        second_regex_samples.append(custom_name)
                    else:
                        second_regex_samples.append(full_name)

        if len(regex_samples) == 0:
            if self.verbose > 0:
                print(f"No samples found with regex '{regex}'"
                      f'{(" after excluding " + str(n_exclude)) if n_exclude > 0 else ""}')
            # Nothing found so return empty array
            return np.array([])

        if create_pairs:
            # Combine sample names into pairs
            regex_sample_pairs = []

            if unique_second_regex:
                if len(second_regex_samples) == 0:
                    if self.verbose > 0:
                        print("No samples found with second regex '" + second_regex + "'")
                    return np.array([])
                
                for sample_1 in regex_samples:
                    for sample_2 in second_regex_samples:
                        if sample_1 != sample_2:
                            regex_sample_pairs.append([sample_1, sample_2])
            else:
                if len(regex_samples) == 1:
                    if self.verbose > 0:
                        print("Cannot form pairs for regex '" + regex + "' as only one sample found")
                    return np.empty((0, 2))
                else:
                    for i, sample_1 in enumerate(regex_samples):
                        for sample_2 in regex_samples[(i + 1):]:
                            regex_sample_pairs.append([sample_1, sample_2])

            return np.array(regex_sample_pairs)
        else:
            return np.array(regex_samples)
        
    def openDefaultSignal(self, bw_idx, chromosome, signal_type = ""):
        raise ValueError("openDefaultSignal not specified")

    def saveBigWig(self, bw_idx = None, custom_sample_name = "", file_name = "", directory = "", 
                   signal_type = "", signals = [], chrom_sizes = {}):
        """
        Save a signal to a bigWig file.

        params: 
            bw_idx:             Index of sample to save to file.
            custom_sample_name: If file name not set, can specify an alternative name for the same.
            file_name:          Name of bigWig file to save signal to.
            directory:          Full directory with folder to save the file to.
            signal_type:        Can be set to specify which type of signal to read.
            signals:            Allows specific signals to be saved to bigWig if specified.
                                If this parameter is set, it must be formatted as a list of signal 
                                arrays/lists per chromosome. e.g. for chromosomes 
                                "chr1", "chr6", "chrX", then signals could be set as
                                [np.array([0.5, 0.8, 0.7, ... , 4.5, 0, 0.1]), # chr1
                                 np.array([0.0, 0.0, 1.5, ... , 0.3, 0, 0.0]) # chr6
                                 np.array([8.2, 3.5, 3.5, ... , 2.7, 4.6, 6.6]])] # chrX
            chrom_sizes:        A dictionary of chromosomes to save signal for, along with their length.
                                e.g. for hg19 {"chr1": 249250621, "chr6": 171115067, "chrX": 155270560}.
        """

        n_signals = len(signals)

        if (bw_idx is None):
            if n_signals == 0:
                raise ValueError("bw_idx must be specified to know which sample to save signal for")
            if custom_sample_name == "":
                # Unknown sample name
                custom_sample_name = "sample"

            # Check to see if user specified specific chromosomes to save
            if len(chrom_sizes) == 0:
                # Default to use chromosomes set initially
                chrom_sizes = self.chrom_sizes

        elif len(chrom_sizes) == 0:
            # Get sample specific chromosome sizes
            bigwig = pyBigWig.open(self.bigwig_paths[bw_idx])
            chrom_sizes = bigwig.chroms()
            chrom_sizes = {c: chrom_sizes[c] for c in self.chromosomes if c in chrom_sizes}

        if custom_sample_name == "":
            sample_name = list(self.sample_names.values())[bw_idx]
        else:
            sample_name = custom_sample_name

        if file_name == "":
            raise ValueError("saveBigWig requires a file name, but none was given")
        
        if directory == "":
            raise ValueError("saveBigWig requires a directory, but none was given")
            
        if not(file_name.endswith(".bw") or file_name.endswith(".bigWig")):
            file_name = file_name + ".bw"
        file_name = os.path.join(directory, file_name)

        if self.verbose > 0:
            print(f'Saving signal for {sample_name} to bigWig "{file_name}"')

        # Create file if doesn't exist
        open(file_name, "w").close()

        # Write signal to bigWig
        bigwig = pyBigWig.open(file_name, "w")

        # Add chromosome sizes to header
        bigwig.addHeader([(chrom, size) for chrom, size in chrom_sizes.items()])
        chromosomes = np.array(list(chrom_sizes.keys()))
        n_chroms = len(chromosomes)

        if n_signals > 0:
            if n_chroms != n_signals:
                raise ValueError(f"{n_signals} signals could not be matched to {n_chroms} chromosomes")

        chroms_added = []

        for chrom_idx in range(n_chroms):
            chrom = chromosomes[chrom_idx]

            if n_signals == 0:
                chrom_signal = self.openDefaultSignal(bw_idx = bw_idx,
                                                      chromosome = chrom, 
                                                      signal_type = signal_type)
                
                if len(chrom_signal) == 0:
                    # Skip the chromosome if signal not found
                    continue
            else:
                try:
                    # Use the signal a user has specified and index the array matching the chromosome
                    chrom_signal = np.array(signals[chrom_idx], dtype = np.float32)
                except:
                    raise ValueError(f"Could not find signal for {chrom} in provided signal list")

            chrom_signal_size = len(chrom_signal)
            
            if chrom_signal_size == 0:
                if self.verbose > 0:
                    print(f"Signal is missing for {sample_name} at {chrom}")
                continue

            elif chrom_signal_size != chrom_sizes[chrom]:
                raise ValueError(f"{chrom} has length {chrom_sizes[chrom]}, but signal has length "
                                 f"{chrom_signal_size}")

            # Remove adjacent duplicated values for saving bigWig intervals
            # e.g. for signal [0.1, 0.1, 0.1, 5.7, 5.7, 5.7, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4, 5.7, ...],
            # extract the values [0.1, 5.7, 3.4, 5.7, ...] and coordinates [0, 3, 6, 12, ...]
            non_consecutive_idxs = np.ones(chrom_signal_size, dtype = bool)
            non_consecutive_idxs[1:] = chrom_signal[1:] != chrom_signal[:-1]
            interval_values = chrom_signal[non_consecutive_idxs]
            interval_coords = np.where(non_consecutive_idxs)[0]
            # Add intervals of normalised signal to new bigWig
            bigwig.addEntries(np.repeat(chrom, len(interval_values)),
                              starts = interval_coords,
                              ends = np.append(interval_coords[1:], len(chrom_signal)),
                              values = interval_values)

            chroms_added.append(chrom)

        # Close to save contents of file
        bigwig.close()

        if (self.verbose > 0) and (len(chroms_added) == 0):
            print(f"Warning: No chromosome signal saved to bigWig {file_name}")

    @staticmethod
    def saveBED(file_name, contents_df, track_name, description, use_score = 0):
        """
        Write contents within a dataframe to BED file.

        params:
            file_name:   Name of BED file to save to.
            contents_df: DataFrame to write as BED.
            track_name:  Custom track name.
            description: Custom description.
            use_score:   Set as 0 for no score, or 1 to use score.
        """

        # Open the file and set to write
        bed_file = open(file_name, "w")

        # Add BED file header
        header = f'track name="{track_name}" '
        header += f'description="{description}" '
        header += f"useScore={use_score}\n"
        bed_file.write(header)
        bed_file.close()

        # Append lines to BED file
        contents_df.to_csv(file_name, sep = "\t", header = False, index = False, mode = "a")

    def createChromSizes(self, chrom_sizes_file = ""):
        """
        Create a tab separated file with chromosome names and their sizes.

        params:
            chrom_sizes_file: File name to save to.
        """

        if not chrom_sizes_file:
            os.makedirs(self.output_directories["bed"], exist_ok = True)
            chrom_sizes_file = os.path.join(self.output_directories["bed"],
                                            f'{self.analysis_name.replace(" ", "_")}.chrom.sizes')
            
        chrom_sizes_df = pd.DataFrame({"chrom": self.chrom_sizes.keys(), 
                                       "size": self.chrom_sizes.values()})
        chrom_sizes_df.to_csv(chrom_sizes_file, header = False, index = False, sep = "\t")
        self.chrom_sizes_file = chrom_sizes_file

    @staticmethod
    def saveBigBed(bed_file, chrom_sizes_file, subprocess_verbose = None):
        """
        Convert a BED file to a bigBed file.

        params:
            chrom_sizes_file:   Path to file of tab separated chromosomes and sizes.
            subprocess_verbose: Verbose of stdout when running subprocess.
        """

        bigbed_file = f"{os.path.splitext(bed_file)[0]}.bb"
        command = f"bedToBigBed {bed_file} {chrom_sizes_file} {bigbed_file}"

        # Run bedToBigBed
        subprocess.run(command, shell = True, stdout = subprocess_verbose)

    def createZoneBED(self, zone_type, create_bigbed = False, chrom_sizes_file = None, 
                      name_postfix = "", desc_postfix = "", file_postfix = ""):
        """
        Save zones to BED file for visualisation in a genome browser.

        params:
            zone_type:        Set as either 'padded' or 'unpadded'.
            create_bigbed:    Set as True to also save a bigBed file.
            chrom_sizes_file: Path to file of tab separated chromosomes and sizes for creating bigBed.
            name_postfix:     Postfix text to add to end of name in BED header.
            desc_postfix:     Postfix text to add to end of description in BED header.
            file_postfix:     Postfix text to add to end of BED file name.
        """

        chrom_col, start_col, end_col = self.region_coords_cols
        sample_names = np.array(list(self.sample_names.values()))

        zone_type = zone_type.lower()
        if zone_type not in ["padded", "unpadded"]:
            raise ValueError(f"Unsupported zone_type {zone_type}.")

        if create_bigbed:
            if chrom_sizes_file is None:
                chrom_sizes_file = os.path.join(self.output_directories["bed"],
                                                f'{self.analysis_name.replace(" ", "_")}.chrom.sizes')

                if not os.path.exists(chrom_sizes_file):
                    self.createChromSizes(chrom_sizes_file)

            else:
                if os.path.exists(chrom_sizes_file):
                    raise FileNotFoundError(f'chrom_sizes_file "{chrom_sizes_file}" does not exist')

            if self.verbose > 0:
                subprocess_verbose = None
            else:
                # Disable any printing
                subprocess_verbose = subprocess.DEVNULL


        if file_postfix:
            file_postfix = file_postfix.replace(" ", "_")

        if self.verbose > 0:
            if create_bigbed:
                print(f"Saving zones to BED and bigBed files")
            else:
                print(f"Saving zones to BED files")

        results_bed_dir = os.path.join(self.output_directories["bed"], "Zones")

        # Create new directory to store BED files
        os.makedirs(results_bed_dir, exist_ok = True)

        # Create empty dictionary to record each zone coordinates per sample
        all_sample_zones = {}
        all_sample_chroms = {}

        for bw_idx in self.sample_ids:
            all_sample_zones[bw_idx] = np.empty((0,2), dtype = np.uint32)
            all_sample_chroms[bw_idx] = {"chroms": [], "n_repeats": []}

        for chrom in self.chromosomes:
            chrom_signal_zones = self.getSampleZones(chrom)[zone_type]
            # chrom_merged_zones = self.getMergedZones(chrom)

            if len(chrom_signal_zones) > 0:
                for bw_idx in self.sample_ids:
                    # Extract the zones for a specific sample for the chromosome
                    sample_zones = chrom_signal_zones[bw_idx]
                    n_sample_zones = len(sample_zones)

                    # Add zones to the dictionary
                    all_sample_zones[bw_idx] = np.concatenate((all_sample_zones[bw_idx], sample_zones))
                    # Record the chromosome name and the number of zones
                    all_sample_chroms[bw_idx]["chroms"].append(chrom)
                    all_sample_chroms[bw_idx]["n_repeats"].append(n_sample_zones)

            elif self.verbose > 0:
                print(f"Warning: No zones found for {chrom}")

        if all(zones.shape[0] == 0 for zones in all_sample_zones.values()):
            raise ValueError("No zones found for any chromosome")

        for bed_idx in self.sample_ids:
            # Only create a BED file if values are present
            write_bed = False
            # Use single sample's name
            sample_name = sample_names[bed_idx]
            # Get zone coordinates for the sample
            all_coords = all_sample_zones[bed_idx]

            if all_coords.shape[0] > 0:
                write_bed = True
                # Get values to create the chromosome column
                bed_chroms = all_sample_chroms[bed_idx]["chroms"]
                n_chrom_repeats = all_sample_chroms[bed_idx]["n_repeats"]

            if write_bed:
                # Set the file name to save to
                bed_file = f"{zone_type}_zones_{sample_name}"
                if file_postfix:
                    bed_file += f"_{file_postfix}.bed"
                else:
                    bed_file += ".bed"
                bed_file_path = os.path.join(results_bed_dir, bed_file)

                # Set the header
                track_name = sample_name

                if name_postfix:
                    track_name += f" {name_postfix}"

                description = f"{zone_type.title()} zones for {sample_name}"

                if desc_postfix:
                    description += f" {name_postfix}"

                # Format the coordinates as a DataFrame
                result_row_names = [f"{sample_name}_{zone_type}_zones_{i}" for 
                                    i in range(1, all_coords.shape[0] + 1)]
                contents_df = pd.DataFrame({chrom_col: np.repeat(bed_chroms, n_chrom_repeats),
                                            start_col: all_coords[:,0],
                                            end_col: all_coords[:,1],
                                            "name": result_row_names})

                # Write results to BED file
                self.saveBED(file_name = bed_file_path, 
                             contents_df = contents_df,
                             track_name = track_name,
                             description = description,
                             use_score = False)

                if create_bigbed:
                    self.saveBigBed(bed_file = bed_file_path,
                                    chrom_sizes_file = chrom_sizes_file,
                                    subprocess_verbose = subprocess_verbose)
                    
            elif self.verbose > 0:
                print(f"Skipping creating BED file for {sample_name} as no {zone_type} zones were found")

    def findRegionsInRange(self, regions, start_coord, end_coord, cap_regions = False):
        """
        Search a set of coordinates to find those within a specified range.

        params:
            regions:     List of arrays of region coordinates.
            start_coord: Region start coordinate (zero indexed).
            end_coord:   Region end coordinate (zero indexed).
            cap_regions: Set as True to reduce the end coordinates to fit within the range. 
                         Leaving as False will leave the coordinates unchanged, but the edges 
                         may exceed the range.
        """

        if start_coord > end_coord:
            raise ValueError("The start coordinate cannot be greater than the end coordinate")

        range_regions = []

        for sample_idx in range(len(regions)):
            # Get regions for the sample
            sample_regions = np.copy(regions[sample_idx])

            # Find regions fully within in coordinates
            inclusive_mask = ((sample_regions[:, 0] < end_coord) &
                              (sample_regions[:, 1] > start_coord))
            # Find regions on edge(s) of coordinates
            edge_mask = (((sample_regions[:, 0] < start_coord) & # Region covers both edges
                          (sample_regions[:, 1] > end_coord)) |
                         ((sample_regions[:, 0] < start_coord) & # Left edge
                          (sample_regions[:, 1] >= end_coord)) |
                         ((sample_regions[:, 0] <= start_coord) & # Right edge
                          (sample_regions[:, 1] > end_coord)))

            # Filter to get regions within range
            sample_regions = sample_regions[inclusive_mask | edge_mask]
            
            if cap_regions and (sample_regions.shape[0] > 0):
                # Cap edge regions to the coordinate range
                sample_regions[0, 0] = max(start_coord, sample_regions[0, 0])
                sample_regions[-1, 1] = min(end_coord, sample_regions[-1, 1])

            # Reassign the regions
            range_regions.append(sample_regions)

        return range_regions

    def plotTracks(self, start_coord, end_coord, chromosome, signals = {}, 
                   peaks = np.array([]), y_intercept = None, y_intercept_label = "",
                   overlay_plots = False, plot_samples = np.array([]), main_title = "", 
                   legend_title = None, custom_colours = np.array([]), peak_colour = "darkgrey", 
                   same_y_axis = True, custom_y = np.array([]), hide_multiple_x_axis = False, 
                   signal_transparency = 0.5, plot_width = 15, plot_height = 0):
        """ 
        Plot signal over a specified region for one or more samples.
        
        params:
            start_coord:          Position of first base pair to plot.
            end_coord:            Position of last base pair to plot.
            chromosome:           Chromosome to plot, e.g. "chr1".
            signals:              Whole chromosome signals can be provided manually by setting 
                                  this as a dictionary of names and signal arrays, or read 
                                  automatically if not set.
            peaks:                List of arrays of coordinates of peaks / regions to plot as bars at 
                                  the top.
            y_intercept:          Set as a number to plot a y-intercept (horizontal line).
            overlay_plots:        Set as True to plot sample signals on top of one another.
            plot_samples:         Array of sample names to plot.
            main_title:           Custom title to display above plot.
            legend_title:         Custom title for legend.
            custom_colours:       List of colours to match to samples.
            peak_colour:          Colour of peak bars (if peaks given).
            same_y_axis:          Whether to plot all samples with the same y-axis range.
            custom_y:             Optionally specify the y-axis range in format [first, last], 
                                  e.g. [0, 200].
            hide_multiple_x_axis: Set as True to only display the x-axis labels for the last sub-plot.
            signal_transparency:  Opacity of the sample's signal.
            plot_width:           Length of the plot.
            plot_height:          Can be used to manually set the height of the plot. 
                                  By default, this is determined by the number of sub-plots.
        """

        if len(signals) == 0:
            custom_sample_names = []

            if len(self.bigwig_paths) == 0:
                raise ValueError("Cannot plot sample signals as bigwig_paths was not set")

            # Get names of all samples
            if isinstance(self.sample_names, dict):
                sample_names = list(self.sample_names.keys())
                custom_sample_names = list(self.sample_names.values())
            elif isinstance(self.sample_names, np.array):
                sample_names = self.sample_names
            elif self.sample_names is None:
                raise ValueError("sample_names was not set")
            else:
                raise ValueError(f"sample_names has unsupported type {type(self.sample_names)}")
            
            if len(plot_samples) > 0:
                if isinstance(plot_samples, str):
                    plot_signals = np.array([plot_samples])
                else:
                    plot_signals = np.array(plot_samples)

                # Check that samples to plot match samples
                if len(np.setdiff1d(plot_samples, sample_names)) > 0:
                    # If avaliable, check custom sample names too
                    if len(custom_sample_names) == 0 or len(np.setdiff1d(custom_sample_names, 
                                                                         sample_names)) > 0:
                        raise Exception("Cannot find samples specified in plot_samples")
            else:
                # Default to plot all samples
                plot_signals = sample_names

            # Convert to indexes to match samples to bigWigs
            plot_bw_idxs = np.array([np.where(np.atleast_1d(sample_names) == s)[0][0] for 
                                     s in plot_signals], dtype = np.uint16)

        elif not isinstance(signals, dict):
            raise ValueError("signals must be set as a dictionary of signal names and arrays, or "
                             "left empty.\n"
                             'For example, signals = {"signal_a": [0,1,2,...,1], "signal_b": '
                             "[1,2,2,...,0]}.")

        else:
            # Derive signal names from signals dictionary
            plot_signals = list(self.signals.keys())
        
        plot_signal_idxs = np.arange(len(plot_signals))

        if not isinstance(self.chrom_sizes, dict):
            raise ValueError("chrom_sizes was not set")

        # Cast coordinates as integers
        start_coord = int(start_coord)
        end_coord = int(end_coord)

        # Error checks for coordinates
        if start_coord > end_coord:
            raise Exception("Start coordinate cannot exceed the end coordinate")
        elif end_coord > self.chrom_sizes[chromosome]:
            raise Exception(f"Cannot create plot as {chromosome} has length "
                            f"{self.chrom_sizes[chromosome]}, but the given end coordinate "
                            f"{end_coord} exceeds this.")

        n_peaks = len(peaks)

        if n_peaks > 0:
            if np.any([isinstance(p, pd.DataFrame) for p in peaks]):
                print("DataFrame peaks not supported")
                return None
            if not isinstance(peaks, list):
                raise ValueError("peaks must be given as a list of arrays")
            elif n_peaks > 1:
                if (not same_y_axis) or (len(peaks) != len(plot_signals)):
                    raise ValueError("Number of peaks and number of samples to plot does not match")

        if isinstance(y_intercept, (int, float)):
            # Plot a y-axis line
            plot_y_intercept = True
        else:
            plot_y_intercept = False

        if len(custom_colours) == len(plot_signals):
            # Create dictionary mapping custom colours to samples to plot
            trace_colour_dict = dict(zip(plot_signals, custom_colours))
        else:
            if self.verbose > 0 and len(custom_colours) > 0:
                print(f"Warning: cannot map {len(custom_colours)} "
                      f"custom colours to {len(plot_signals)} samples")
            # Create list of hex colours to represent all samples
            trace_colours = np.array(list(color_palette("bright", len(sample_names)).as_hex()))
            # Convert to dictionary mapping samples to plot with their colour
            trace_colour_dict = dict(zip(plot_signals, [trace_colours[i] for i in plot_signal_idxs]))

        coords_diff = end_coord - start_coord

        # If signals not provided, read signal from file
        if len(signals) == 0:
            signals = np.zeros((len(plot_signals), coords_diff))

            if (len(plot_signals) == 1) or (self.n_cores == 1):
                bw_idx = plot_bw_idxs[0]
                signals[0] = self.extractChunkSignal(bw_idx = bw_idx,
                                                     bigwig_file = self.bigwig_paths[bw_idx], 
                                                     chromosome = chromosome,
                                                     start_idx = start_coord,
                                                     end_idx = end_coord,
                                                     pad_end = True)

            else:
                with ProcessPoolExecutor(self.n_cores) as executor:
                    read_signal_processes = [executor.submit(self.extractChunkSignal,
                                                             bw_idx = bw_idx,
                                                             bigwig_file = self.bigwig_paths[bw_idx], 
                                                             chromosome = chromosome,
                                                             start_idx = start_coord,
                                                             end_idx = end_coord,
                                                             pad_end = True,
                                                             return_file = True) for bw_idx in plot_bw_idxs]
                
                for process in as_completed(read_signal_processes):
                    # Get the signal and the file as an identifier
                    bigwig_file, signal = process.result()
                    # Convert file into signal index
                    bw_idx = np.where(self.bigwig_paths == bigwig_file)[0][0]
                    row_idx = np.where(plot_bw_idxs == bw_idx)[0][0]
                    signals[row_idx] = signal
        else:
            # Convert to numpy array
            signals = np.vstack(signals.values())
            signal_length = signals.shape[0]

            if signal_length < coords_diff:
                raise ValueError(f"Signals are smaller than range {chromosome}:{start_coord}-{end_coord}")
            elif signal_length != end_coord - start_coord:
                # If full signal given, sub-set
                signals = signals[:, start_coord:end_coord]

        plot_peaks = False

        if len(peaks) > 0:
            # Filter to find those within the coordinates and cap their boundaries
            peaks = self.findRegionsInRange(regions = peaks,
                                            start_coord = start_coord,
                                            end_coord = end_coord,
                                            cap_regions = True)
            
            for i in range(len(peaks)):
                peak_coords = peaks[i]
                if not plot_peaks:
                    if np.sum(peak_coords) > 0:
                        # Found one or more valid peaks to plot
                        plot_peaks = True
            
        border_pad = 0
        use_global_y = True

        if len(custom_y) == 2:
            min_y, max_y = custom_y
        elif overlay_plots or same_y_axis:
            # Calculate the global lowest and highest signal intensity
            min_y = np.nanmin(signals)
            max_y = np.nanmax(signals)
            # Round to nearest 5
            min_y = round(math.floor(min_y / 5)) * 5
            max_y = round(math.ceil(max_y / 5)) * 5
            # Add 1/10th whitespace to top/bottom of plot to prevent signal and boundary overlapping
            if min_y < 0:
                border_pad = np.abs(min_y / 10)
            if max_y > 0:
                border_pad = max((max_y / 10), border_pad)

            if max_y == 0 and min_y == 0:
                print("Warning: All signals were zero")
        else:
            # Each subplot has its own y-axis height
            use_global_y = False

        if use_global_y and plot_peaks:
            peak_boundary_height = max_y / 20

        # Set plot dimensions
        n_plots_above = 0

        if overlay_plots:
            # Plot single graph
            n_subplots = 1 + n_plots_above
            if plot_height <= 0:
                plot_height = 3 + (n_plots_above * 2.5)
            height_rations = (([1] * n_plots_above) + [2.5])
        else:
            # Make plot longer as including multiple subplots
            n_subplots = len(signals) + (len(signals) * n_plots_above)
            if plot_height <= 0:
                plot_height = (3 * len(signals)) + (n_plots_above * 2.5)
            height_rations = (([1] * n_plots_above) + [2.5]) * len(signals)

        # Create subplots
        fig, ax = plt.subplots(nrows = n_subplots, ncols = 1, figsize = (plot_width, plot_height),
                               gridspec_kw = {'height_ratios': height_rations})

        # Set the x-axis coordinates shared for all tracks          
        x = np.arange(start_coord, end_coord)
        last_row_idx = len(signals) - 1

        # Iterate through samples
        for signal_name, row_idx in zip(plot_signals, plot_signal_idxs):
            row_signal = signals[row_idx]
            row_min_y = round(math.floor(min(row_signal) / 5)) * 5
            row_max_y = round(math.ceil(max(row_signal) / 5)) * 5

            # Set the axis to plot the track on
            if n_subplots == 1:
                plot_ax = ax
            else:
                if overlay_plots:
                    plot_ax = ax[n_plots_above]
                else:
                    plot_ax = ax[(n_plots_above * (row_idx + 1)) + row_idx]

            # Plot the y axis coordinates for the specific track
            plot_ax.plot(x,
                         row_signal,
                         color = trace_colour_dict[signal_name],
                         alpha = signal_transparency)

            try:
                if row_min_y < 0:
                    # Fill area above line
                    plot_ax.fill_between(x,
                                         row_signal,
                                         where = row_signal < 0,
                                         color = trace_colour_dict[signal_name],
                                         label = signal_name,
                                         alpha = signal_transparency)
                if row_max_y > 0:
                    # Fill area below line
                    plot_ax.fill_between(x,
                                         row_signal,
                                         where = row_signal > 0,
                                         color = trace_colour_dict[signal_name],
                                         label = signal_name,
                                         alpha = signal_transparency)
                    
            except Exception as e:
                print(f"Error: Could not plot filled area under curve for sample {signal_name}")
                print(e)

            if (len(custom_y) != 2) and (not same_y_axis) and (not overlay_plots):
                # Calculate local lowest and highest signal intensity
                min_y = row_min_y
                max_y = row_max_y

                if (self.verbose > 0) and (max_y == 0) and (min_y == 0):
                    print(f"Warning: the signal for {signal_name} was zero for the region")
                elif min_y > 0:
                    min_y = 0
                elif max_y < 0:
                    max_y = 0

                # Add 1/10th whitespace to top of plot to prevent signal and boundary overlapping
                if min_y < 0:
                    border_pad = np.abs(min_y / 10)
                if max_y > 0:
                    border_pad = max((max_y / 10), border_pad)

            if plot_peaks:
                # Determine which peaks to plot (if any)
                peak_coords = np.empty((0, 2))

                if overlay_plots:
                    if row_idx == 0:
                        # Only plot one lot of peaks if all signal on same axis
                        peak_coords = peaks[0]
                else:
                    if len(peaks) == 1:
                        # Plot the same peaks for each sample
                        peak_coords = peaks[0]
                    else:
                        # Plot peaks specific to the sample
                        peak_coords = peaks[row_idx]

                if not use_global_y:
                    # Use individual y-axis height to set the boundary
                    peak_boundary_height = max_y / 20

                if peak_coords.shape[0] > 0:
                    # Iterate over peaks for sample
                    for coords in peak_coords:
                        peak_x = coords[0]
                        peak_y = coords[1]

                        # Plot peak boundary at the top of the plot
                        plot_ax.add_patch(patches.Rectangle((peak_x, # lower-left x-coord
                                                             max_y - peak_boundary_height + border_pad), # lower-left y-coord
                                                             peak_y - peak_x, # Width
                                                             peak_boundary_height + border_pad, # Height
                                                             color = peak_colour,
                                                             linewidth = 0))

            if (overlay_plots and (row_idx == 0)) or (not overlay_plots):
                # Set the axis limits
                plot_ax.set_ylim([min_y, max_y + border_pad])
                plot_ax.set_xlim([start_coord, end_coord])
                # Set y-axis title
                plot_ax.set_ylabel("Signal")

            if hide_multiple_x_axis and row_idx != last_row_idx:
                # Hide x-axis labels for each sub-plot except the last
                plot_ax.get_xaxis().set_visible(False)

            if plot_y_intercept:
                plot_ax.axhline(y = y_intercept, color = "black", linestyle = "--", 
                                label = y_intercept_label)

            plot_ax.legend(title = legend_title)
            plot_ax.ticklabel_format(style = "plain", useOffset = False)

        if len(main_title) > 0:
            # Set a main title above the subplots
            fig.suptitle(main_title, fontweight = "bold", y = 1.01)

        # Disable scientific notation and display plot
        plt.ticklabel_format(style = "plain", useOffset = False)
        plt.show()