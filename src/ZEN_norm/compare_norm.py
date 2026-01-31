import numpy as np
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap
from concurrent.futures import ProcessPoolExecutor, as_completed
from .chrom_analysis import ChromAnalysisCore

##########################################
# Main class for normalisation comparision
##########################################

class CompareNorm(ChromAnalysisCore):
    __slots__ = ("bigwig_df", "regions_df", "coords_df", "min_peak_score", "min_consensus", "border_pad", 
                 "norm_methods", "chrom_maxs")

    def __init__(self, bigwig_df, regions_df = None, coords_df = None, min_peak_score = None, 
                 min_consensus = 1, border_pad = 0, n_cores = 1, analysis_name = "Analysis", 
                 verbose = 1):
        """
        Class to create plots for comparing normalisation method performance.
        
        params:
            bigwig_df:      How input bigWigs are mapped to the samples they belong to and the normalisation 
                            method that was applied. Settable as either a file path to a CSV or as a 
                            DataFrame with columns "norm", "sample" and "bigwig", where each row contains 
                            the name of the normalisation method, sample name and bigWig path.
            regions_df:     To evaluate normalisation performance over regions such as peaks, zones or 
                            custom coordinates, set as a DataFrame with columns "sample" and "regions", 
                            where regions are given by paths to BED or narrowPeak files.
            coords_df:      To evaluate normalisation performance over custom coordinates, set as a 
                            DataFrame containing "chrom", "start" and "end". Must be given if regions_df 
                            is not specified.
            min_peak_score: Threshold between 0 to 1 to filter peaks by if LanceOtron peaks are given 
                            in regions_df.
            min_consensus:  Minimum number of samples that a region is present in for it to be included 
                            when evaluating each normalisation method. By default this is 1 so that 
                            all regions are included.
            border_pad:     How close region coordinates have to be to be combined, 
                            e.g. [1,000, 3,000] and [4,000, 10,000] with border_pad = 1,000 would be merged.
            n_cores:        The number of cores / CPUs to use if using multiprocessing.
            analysis_name:  Custom name of folder to save results to. By default this will be set to 
                            "Analysis".
            verbose:        Set as an integer greater than 0 to display progress messages.
        """

        # Initialise the parent class
        super().__init__(n_cores = n_cores,
                         analysis_name = analysis_name,
                         verbose = verbose)

        self.output_directories["plots"] = os.path.join(self.output_directories["results"], "Plots")
        self.output_directories["output_stats"] = os.path.join(self.output_directories["results"], "Stats")

        self.setMinPeakScore(min_peak_score)
        self.setMinConsensus(min_consensus)
        self.setBorderPad(border_pad)
        self.setBigWigs(bigwig_df = bigwig_df)
        self.setCoords(regions_df = regions_df, coords_df = coords_df)
        self.setChromAttributes()
        
    def __str__(self):
        """ Format parameters in user readable way when print is called on an object """

        if self.verbose <= 0:
            verbose_msg = "(silent)"
        elif self.verbose == 1:
            verbose_msg = "(active)"
        else:
            verbose_msg = "(debugging mode)"

        min_consensus_msg = f"regions must be present in at least {self.min_consensus} sample"

        if self.min_consensus != 1:
            min_consensus_msg += "s"

        message = (f'{self.__class__.__name__} object for "{self.analysis_name}"\n'
                   f"   * Number of samples: {len(self.sample_names)}\n"
                   f"   * Sample names: {', '.join(self.sample_names)}\n"
                   f"   * Normalisation methods: {', '.join(self.norm_methods)}\n"
                   f"   * LanceOtron minimum peak score: {self.min_peak_score}\n"
                   f"   * Minimum consensus: {min_consensus_msg}\n"
                   f"   * Resources: {self.n_cores} cores\n"
                   f"   * Verbose: {self.verbose} {verbose_msg}")

        return message

    def setNCores(self, n_cores):
        try:
            self.n_cores = max(int(n_cores), 1)
        except:
            raise ValueError("n_cores must be a non-negative integer as this is the number of "
                             "CPUs to use for parallelisation.")

    def setMinPeakScore(self, min_peak_score):
        if min_peak_score is not None:
            error_message = "Minimum peak score must be set as a value between 0 - 1, "
            error_message += "or set as None to disable LanceOtron peak score filtering"

            try:
                min_peak_score = float(min_peak_score)

                if (min_peak_score < 0) or (min_peak_score > 1):
                    raise ValueError(error_message)
                
            except:
                raise ValueError(error_message)
            
        self.min_peak_score = min_peak_score

    def setMinConsensus(self, min_consensus):
        error_message = "Minimum number of samples for consensus regions, min_consensus, must be set as "
        error_message += "a positive integer"
        try:
            min_consensus = int(min_consensus)

            if min_consensus > 0:
                self.min_consensus = min_consensus
            else:
                raise ValueError(error_message)

        except:
            raise ValueError(error_message)

    def setBorderPad(self, border_pad):
        error_message = "Border pad must be set as zero for no padding, or a positive integer"
        try:
            border_pad = int(border_pad)

            if border_pad >= 0:
                self.border_pad = border_pad
            else:
                raise ValueError(error_message)

        except:
            raise ValueError(error_message)

    def setBigWigs(self, bigwig_df):
        """
        Set a DataFrame containing names of normalisation methods, sample names and paths to bigWigs.

        params:
            bigwig_df: DataFrame with columns 'norm', 'sample' and 'bigwig'.
        """

        if isinstance(bigwig_df, str):
            if os.exists(bigwig_df):
                # Read from file
                bigwig_df = pd.read_csv(bigwig_df)
            else:
                raise FileNotFoundError(f'Cannot open "{bigwig_df}" as the file was not found')
                
        elif not isinstance(bigwig_df, pd.DataFrame):
            raise ValueError(f"bigwig_df must be set as either a DataFrame or a path to a csv")
            
        if len(bigwig_df) == 0:
            raise ValueError("bigwig_df is empty")

        columns = []
            
        for col in bigwig_df.columns:
            # Sanitise the column names
            col = col.strip().lower()

            if col.endswith("s"):
                # Remove 's' from end
                col = col[:-1]
            if col.startswith("norm"):
                    # Set variations of 'normalisation' as 'norm'
                    col = "norm"
            if col == "bw":
                # Replace shorthand for bigwig
                col = "bigwig"

            columns.append(col)

        # Check if any expected columns are missing
        expected_columns = ["norm", "sample", "bigwig"]
        missing_columns = np.setdiff1d(columns, expected_columns)

        if len(missing_columns) > 0:
            raise ValueError(f"bigwig_df must contain columns {', '.join(expected_columns)}")

        # Set the sanitised columns
        bigwig_df.columns = columns
        # Keep only relevant columns and ensure each row contains strings
        bigwig_df = bigwig_df[expected_columns].astype(str)

        sample_names = bigwig_df["sample"].unique()

        if len(sample_names) < 2:
            raise ValueError("Only one sample found, but at least two are required")

        # Check that each file is a bigWig
        bigwig_files = bigwig_df["bigwig"].to_numpy().astype(str)
        n_invalid = np.sum(~(np.char.endswith(bigwig_files, ".bw") | 
                             np.char.endswith(bigwig_files, ".bigWig")))

        if n_invalid > 0:
            raise ValueError(f"{n_invalid} file path{'s' if n_invalid != 1 else ''} were missing "
                             f"a bigwig file extension")

        self.bigwig_df = bigwig_df
        self.norm_methods = bigwig_df["norm"].unique()
        self.sample_names = self.naturalSort(sample_names)

    def setCoords(self, regions_df = None, coords_df = None):
        """
        Set the coordinates DataFrame from either a given DataFrame or region BED files. 
        Any overlapping coordinates are merged.

        params:
            regions_df: DataFrame containing 'sample' and 'regions', where region coordinates are
                        given by paths to BED or narrowPeak files.
            coords_df:  Non-sample specific DataFrame containing 'chrom', 'start' and 'end'.
        """

        if isinstance(regions_df, pd.DataFrame):
            # Check expected columns exist
            columns = []
                    
            for col in regions_df.columns:
                # Sanitise the column names
                col = col.strip().lower()
                
                if col.startswith("region") or col.startswith("zone") or col.startswith("coord"):
                    # Rename column for consistency
                    col = "regions"

                columns.append(col)

            # Check if any expected columns are missing
            expected_columns = ["sample", "regions"]
            missing_columns = np.setdiff1d(columns, expected_columns)

            if len(missing_columns) > 0:
                raise ValueError(f"regions_df must contain columns {', '.join(expected_columns)}")

            # Keep only relevant columns
            regions_df.columns = columns
            regions_df = regions_df[expected_columns]

            missing_samples = np.setdiff1d(self.sample_names, regions_df["sample"])
            n_missing = len(missing_samples)

            if n_missing > 0:
                raise ValueError(f"{n_missing} samples were missing files in the regions_df")
            
            regions_paths = regions_df["regions"].to_list()

            for file in regions_paths:
                if (not file.endswith(".bed")) and (not file.endswith(".narrowPeak")):
                    raise ValueError(f'Regions file "{file}" is not a BED or narrowPeak file')
                elif not os.path.exists(file):
                    raise FileNotFoundError(f'Regions file "{file}" does not exist')
                
            # Set coordinates as merged coordinates across all regions
            self.regions_df = regions_df
            self.coords_df = self.mergeRegions(regions_paths, min_consensus = self.min_consensus)

        elif isinstance(coords_df, pd.DataFrame):
            if coords_df.empty:
                raise ValueError("coords_df is empty")
            
            # Check chromosome coordinate columns exist
            columns = []
                
            for col in coords_df.columns:
                # Sanitise the column names
                col = col.strip().lower()
            
                if col.startswith("chrom"):
                    col = "chrom"

                columns.append(col)

            # Check if any expected columns are missing
            expected_columns = ["chrom", "start", "end"]
            missing_columns = np.setdiff1d(columns, expected_columns)

            if len(missing_columns) > 0:
                raise ValueError(f"coords_df must contain columns {', '.join(expected_columns)}")

            # Keep only relevant columns
            coords_df.columns = columns
            coords_df = coords_df[expected_columns]

            chrom_coords = {}

            # Combine any overlapping coordinates
            for chrom, chrom_rows in coords_df.groupby("chrom"):
                chrom_rows = chrom_rows[["start", "end"]].to_numpy()
                chrom_coords[chrom] = self.mergeOverlapCoords(chrom_rows)

            # Convert into dataframe
            self.coords_df = pd.DataFrame({"chrom": np.hstack([[chrom] * len(chrom_coords[chrom]) for 
                                                               chrom in chrom_coords]),
                                           "start": np.hstack([chrom_coords[chrom][:,0] for 
                                                               chrom in chrom_coords]),
                                           "end": np.hstack([chrom_coords[chrom][:,1] for 
                                                             chrom in chrom_coords])})
            self.regions_df = None
            self.min_consensus = 1

        else:
            raise ValueError("regions_df or coords_df must be set as a DataFrame of region files or coordinates")

    def setChromAttributes(self):
        """
        Set chromosomes and their sizes.
        """

        self.chromosomes = self.orderChroms(chromosomes = np.unique(self.coords_df["chrom"]))

        # Open each bigWig file specified in the dataframe
        bigwig_files = self.bigwig_df["bigwig"].to_numpy().astype(str)
        bigwigs, _ = self.openFiles(file_paths = bigwig_files, verbose = self.verbose)

        # Get the size of each chromosome
        self.chrom_sizes = self.calculateChromSizes(bigwigs = bigwigs, chromosomes = self.chromosomes)

        # Check that coordinates match the chromosome sizes
        for chrom in self.chromosomes:
            chrom_coords = np.array(self.coords_df[self.coords_df["chrom"] == chrom][["start"]])

            if np.any(chrom_coords > self.chrom_sizes[chrom]):
                raise ValueError(f"Coordinates found for {chrom} in coords_df that exceeded the "
                                 f"chromsome size of {self.chrom_sizes[chrom]}")

        # Set the maximum value per chromosome for each sample
        chrom_maxs = {}

        for bw_idx in range(len(bigwigs)):
            bigwig = bigwigs[bw_idx]
            chrom_maxs[bw_idx] = {}

            for chrom in self.chromosomes:
                if bigwig.chroms(chrom) is not None:
                    chrom_maxs[bw_idx][str(chrom)] = bigwig.stats(chrom, type = "max")[0]

        self.chrom_maxs = chrom_maxs

    def getSampleNames(self):
        return self.sample_names

    def regexFindSamples(self, regex, second_regex = "", ignore_case = False, 
                         exclude_samples = np.array([])):
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
            exclude_samples:     Optionally set list of sample names to ignore.
        """

        if ignore_case:
            # Make case-insensitive
            re_case = re.IGNORECASE
        else:
            re_case = 0

        # Get sample names
        sample_names = self.sample_names
        exclude_samples = np.array(exclude_samples)
        n_exclude = 0

        if len(exclude_samples) > 0:
            # Find indexes of the sample names to ignore
            exclude_idxs = np.where(np.isin(sample_names, exclude_samples))[0]
            n_exclude = len(exclude_idxs)

            if n_exclude > 0:
                # Remove names of samples
                sample_names = np.delete(sample_names, exclude_idxs)

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
        for sample in sample_names:
            check = re.search(regex, sample, re_case)

            if check:
                regex_samples.append(sample)

            elif unique_second_regex:
                check = re.search(second_regex, sample, re_case)
                if check:
                    second_regex_samples.append(sample)

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
        
    @staticmethod
    def createBigWigFileCSV(sample_names, bigwig_directory, norm_methods, file_extension = ".bw", 
                            csv_file = ""):
        """
        Create a DataFrame mapping normalised bigWigs with normalisation methods and sample names. 
        bigWig files are assumed to be in subfolders within the given directory and named with convention 
        [bigwig_directory]/[normalisation method]/[sample]_[normalisation method][file_extension]. 
        The exception is for non-normalised (raw) bigWigs, which are named with convention 
        [bigwig_directory]/No_Normalisation/[sample][file_extension].

        params:
            sample_names:     List or array of sample names to gather files for.
            bigwig_directory: Directory containing subfolders with bigWigs per normalisation method.
            norm_methods:     List or array of normalisation methods to collect file paths for.
            file_extension:   Either '.bw' or '.bigWig'.
            csv_file:         Set as a file name to save the DataFrame to a CSV.

        returns:
            bigwig_df: DataFrame with columns 'norm' 'sample' and 'bigwig'.
        """

        if len(sample_names) == 0:
            raise ValueError("No samples given")
            
        if not os.path.isdir(bigwig_directory):
            raise ValueError(f'Invalid path "{bigwig_directory}"')
            
        if not file_extension.startswith("."):
            file_extension = "." + file_extension
        if file_extension not in [".bw", ".bigWig"]:
            raise ValueError('File extension must be either ".bw" or ".bigWig"')

        save_to_csv = False

        if csv_file != "":
            save_to_csv = True

            if not csv_file.endswith(".csv"):
                csv_file = csv_file + ".csv"

        bigwig_files = []

        for norm in norm_methods:
            norm = str(norm)

            if norm.lower() in ["raw", "no_normalisation", "none"]:
                non_normalised = True
            else:
                non_normalised = False

            for sample in sample_names:
                if non_normalised:
                    # Set path for non-normalised (raw) bigWigs
                    bigwig_files.append(f"{bigwig_directory}/No_Normalisation/{sample}{file_extension}")
                else:
                    # Set path for normalised bigWigs
                    bigwig_files.append(f"{bigwig_directory}/{norm}/{sample}_{norm}{file_extension}")

        # Create a CSV with rows per normalisation method, sample and bigWig path
        bigwig_df = pd.DataFrame({"norm": np.array([[m] * len(sample_names) for m in norm_methods]).flatten(), 
                                  "sample": np.tile(sample_names, len(norm_methods)),
                                  "bigwig": bigwig_files})
        
        if save_to_csv:
            bigwig_df.to_csv(csv_file, header = True, index = False)

        return bigwig_df

    @staticmethod
    def createRegionsFileCSV(sample_names, regions_directory, file_prefix = "", file_postfix = "", 
                             file_extension = ".bed", csv_file = ""):
        """
        Create a DataFrame mapping sample names to files with region coordinates (either BED or narrowPeak). 
        Region files are assumed to be in the same directory, and named with convention 
        [regions_directory]/[file_prefix][sample][file_postfix][file_extension].

        params:
            sample_names:       List or array of sample names to gather files for.
            regions_directory:  Directory containing BED or narrowPeak files.
            file_prefix:        Beginning part of the name of each file.
            file_postfix:       End part of the name of each file.
            file_extension:     Type of file region coordinates are saved as. Either '.bed' or '.narrowPeak'.
            csv_file:           Set as a file name to save the DataFrame to a CSV.

        returns:
            regions_df: DataFrame with columns 'sample' and 'regions'.
        """
        
        region_files = []

        if len(sample_names) == 0:
            raise ValueError("No samples given")
        
        if not os.path.isdir(regions_directory):
            raise ValueError(f'Invalid path "{regions_directory}"')
        
        if not file_extension.startswith("."):
            file_extension = "." + file_extension
        if file_extension not in [".bed", ".narrowPeak"]:
            raise ValueError('File extension must be either ".bed" or ".narrowPeak"')

        save_to_csv = False

        if csv_file != "":
            save_to_csv = True

            if not csv_file.endswith(".csv"):
                csv_file = csv_file + ".csv"

        for sample in sample_names:
            # Create list of BED files
            sample_regions_file = os.path.join(regions_directory, 
                                               f"{file_prefix}{sample}{file_postfix}{file_extension}")
            region_files.append(sample_regions_file)

            if not os.path.exists(sample_regions_file):
                print(f"Warning: File {sample_regions_file} does not exist")

        regions_df = pd.DataFrame({"sample": sample_names,
                                   "regions": region_files})

        if save_to_csv:
            regions_df.to_csv(csv_file, header = True, index = False)

        return regions_df

    @staticmethod
    def mergeOverlapCoords(chrom_coords, border_pad = 0, min_consensus = 1, coord_ids = []):
        """ 
        Find chromosome coordinates that overlap and combine them. 
        e.g. if chrom_coords = [[100, 200], [300, 500], [700, 1000], [50, 200], [300, 600]], 
        coord_ids = [1, 1, 1, 2, 2], border_pad = 50 and min_consensus = 2, then the merged coordinates 
        would be [[0, 250], [250, 650]].
        
        params:
            chrom_coords:  2D array containing pairs of start and end coordinates for a chromosome.
            border_pad:    How close coordinates have to be to be combined, e.g. [1,000, 3,000] and 
                           [4,000, 10,000] with border_pad = 1,000 would be merged.
            min_consensus: Minimum number of samples that a peak is present in for it to be included 
                           when evaluating each normalisation method.
            coord_ids:     If merging coordinates found in more than one sample, this must be set as 
                           a list of IDs to distinguish which coordinate belongs to which sample. 
        """

        if len(chrom_coords) == 0:
            return np.empty((0, 2), dtype = np.uint32)

        chrom_coords = np.array(chrom_coords).astype(np.uint32).reshape(-1, 2)
        min_consensus = int(min_consensus)

        # Merge all coordinates
        if min_consensus <= 1:
            # Remove duplicates
            chrom_coords = np.unique(chrom_coords, axis = 0)

            if len(chrom_coords) == 1:
                return chrom_coords

            # Sort by first coordinate
            chrom_coords = chrom_coords[chrom_coords[:, 0].argsort()]

            merged_coords = []
            current_start = chrom_coords[0][0]
            current_end = chrom_coords[0][1]

            for start, end in chrom_coords[1:]:
                # Test if new coordinates falls within the range of the current ones
                if start <= current_end + border_pad:
                    current_end = max(current_end, end)
                else:
                    # Record the old coordinate and update the current one
                    merged_coords.append([current_start, current_end])
                    current_start = start
                    current_end = end

            # Add the last coordinate
            merged_coords.append([current_start, current_end])

        # Merge coordinates for those found in a minimum number of samples
        else:
            coord_ids = np.array(coord_ids)
            n_coord_ids = len(coord_ids)
            n_chrom_coords = len(chrom_coords)
            n_unique_coord_ids = len(np.unique(coord_ids))

            if n_coord_ids != len(chrom_coords):
                raise ValueError(f"Merging coordinates that overlap at least min_consensus "
                                 f"{min_consensus} samples requires the number of coordinate IDs "
                                 f"({n_coord_ids}) to match the number of coordinates "
                                 f"({n_chrom_coords})")
            if n_unique_coord_ids < min_consensus:
                raise ValueError(f"Minimum number of samples to overlap, min_consensus ({min_consensus}), "
                                 f"exceeded number of unique coordinates IDs ({n_unique_coord_ids})")

            # Sort by first coordinate
            sort_idxs = chrom_coords[:, 0].argsort()
            chrom_coords = chrom_coords[sort_idxs]
            coord_ids = np.array(coord_ids, dtype = np.uint16)[sort_idxs]

            # Initialise group end as the first coordinate's end
            current_end = chrom_coords[0][1]
            # Check coordinates in groups where they overlap
            group_ids = set([coord_ids[0]])
            group_coords = [chrom_coords[0]]
            merged_coords = []

            for coord_id, coords in zip(coord_ids[1:], chrom_coords[1:]):
                start = coords[0]
                end = coords[1]

                # Test if new coordinates falls within the range of the current ones
                if start <= current_end + border_pad:
                    group_ids.add(coord_id)
                    group_coords.append(coords)
                    # Update the end of the group
                    current_end = max(current_end, end)

                else:
                    if len(group_ids) >= min_consensus:
                        # Merge coordinates only if a minimum number of samples is reached
                        group_merged = CompareNorm.mergeOverlapCoords(group_coords, 
                                                                      border_pad = border_pad,
                                                                      min_consensus = 1)
                        merged_coords.extend(group_merged.tolist())
                
                    # No overlap, so start new group
                    group_ids = set([coord_id])
                    group_coords = [coords]
                    current_end = end

            # Check last group
            if len(group_ids) >= min_consensus:
                group_merged = CompareNorm.mergeOverlapCoords(group_coords,
                                                              border_pad = border_pad,
                                                              min_consensus = 1)
                merged_coords.extend(group_merged.tolist())

        # Convert output to array
        merged_coords = np.array(merged_coords, dtype = np.uint32).reshape(-1, 2)

        return merged_coords

    def mergeRegions(self, regions_paths, min_consensus = 1, border_pad = 0):
        """
        Given a list of region coordinate files, combine the coordinates within them.

        params:
            regions_paths: List of file paths to region coordinates in BED or narrowPeak format.
            min_consensus: Minimum number of samples that a region is present in for it to be included 
                           when evaluating each normalisation method.
            border_pad:    How close coordinates have to be to be combined.
                           e.g. [1,000, 3,000] and [4,000, 10,000] with border_pad = 1,000 would be 
                           merged.
        """

        chrom_coords = {}
        coord_file_idxs = {}

        # Get coordinates for each chromosome across all peaks
        for file_idx, file in enumerate(regions_paths):
            if os.path.exists(file):
                if file.endswith(".bed"):
                    # Try read LanceOtron peaks
                    region_coords = pd.read_csv(file, sep = "\t", header = 0)

                    if len(region_coords.columns) == 1:
                        # Header found in BED file, so try again
                        region_coords = pd.read_csv(file, sep = "\t", header = None, skiprows = 1)

                    elif (self.min_peak_score is not None) and ("overall_peak_score" in region_coords.columns):
                        # Filter peaks by peak score
                        region_coords = region_coords[region_coords["overall_peak_score"] >= self.min_peak_score]

                elif file.endswith(".narrowPeak"):
                    region_coords = pd.read_csv(file, sep = "\t", header = None)

                else:
                    raise ValueError(f'Unsupported file extension for "{file}"')

                # Keep only coordinates
                region_coords = region_coords.iloc[:, :3]
                region_coords.columns = ["chrom", "start", "end"]

            else:
                raise FileNotFoundError(f'Regions file "{file}" does not exist')

            for chrom, chrom_rows in region_coords.groupby("chrom"):
                chrom_rows = chrom_rows[["start", "end"]].to_numpy()

                if chrom in chrom_coords:
                    chrom_coords[chrom].append(chrom_rows)
                else:
                    chrom_coords[chrom] = [chrom_rows]

                if min_consensus > 1:
                    # Record file IDs to keep coordinates that overlap with a minimum number of samples
                    if chrom in coord_file_idxs:
                        coord_file_idxs[chrom].extend([file_idx] * len(chrom_rows))
                    else:
                        coord_file_idxs[chrom] = [file_idx] * len(chrom_rows)

        if len(chrom_coords) == 0:
            raise ValueError("No region coordinates found")

        merged_chrom_coords = {}
        coord_ids = []
        merge_chrom = True

        for chrom in chrom_coords:
            if min_consensus > 1:
                merge_chrom = False

                # Check that regions were found from the minimum number of samples
                if chrom in coord_file_idxs:
                    coord_ids = coord_file_idxs[chrom]

                    if len(np.unique(coord_ids)) >= min_consensus:
                        # Allow merging as enough samples found
                        merge_chrom = True

            if merge_chrom:
                # Combine overlapping coordinates
                merged_coords = self.mergeOverlapCoords(chrom_coords = np.vstack(chrom_coords[chrom]), 
                                                        min_consensus = min_consensus,
                                                        coord_ids = coord_ids,
                                                        border_pad = border_pad)
                if len(merged_coords) > 0:
                    merged_chrom_coords[chrom] = merged_coords

        del chrom_coords

        if len(merged_chrom_coords) == 0:
            if self.verbose > 0:
                print(f"Warning: no regions found after merging")

            coords_df = pd.DataFrame({"chrom": [], "start": [], "end": []})

        else:
            # Convert coordinates into dataframe
            coords_df = pd.DataFrame({"chrom": np.hstack([[chrom] * len(merged_chrom_coords[chrom]) for 
                                                          chrom in merged_chrom_coords]),
                                      "start": np.hstack([merged_chrom_coords[chrom][:,0] for 
                                                          chrom in merged_chrom_coords]),
                                      "end": np.hstack([merged_chrom_coords[chrom][:,1] for 
                                                        chrom in merged_chrom_coords])})

        return coords_df

    def signalReader(self, bw_idx, norm_method, chromosome, start_idx = 0, end_idx = -1, pad_end = False,
                     dtype = np.float32, verbose = 0):
        """
        Read bigWig signal into a numpy array.

        params:
            bw_idx:      Index of sample to read signal for.
            norm_method: Name of the normalisation method in bigwig_df to read signal for.
            chromosome:  Chromosome to read signal for.
            start_idx:   Start base-pair position (zero indexed).
            end_idx:     End base-pair position (zero indexed).
            pad_end:     Whether to add zeros to end to ensure array size is consistent across 
                         samples for a chromosome
            dtype:       Data type to store the signal as.
            verbose:     Set as a number > 0 to print progress.
        """

        if norm_method not in self.norm_methods:
            raise ValueError(f'Unknown normalisation method "{norm_method}"')

        sample_name = self.sample_names[bw_idx]
        bigwig_file = self.bigwig_df.loc[(self.bigwig_df["norm"] == norm_method) &
                                         (self.bigwig_df["sample"] == sample_name), "bigwig"].squeeze()
        
        # Open the bigWig for the sample
        return self.extractChunkSignal(bigwig_file = bigwig_file, 
                                       sample_name = sample_name,
                                       chromosome = chromosome,
                                       start_idx = start_idx,
                                       end_idx = end_idx,
                                       pad_end = pad_end,
                                       dtype = dtype,
                                       verbose = verbose)

    def sampleCoordCounts(self, bw_idx, chromosome, coords, norm_method):
        """
        Given a set of coordinates, calculate the sum of the signal across each one.

        params:
            bw_idx:      Sample ID to read signal for.
            chromosome:  Name of chromosome to read signal for. 
            coords:      Array of start and end coordinate pairs.
            norm_method: Normalisation method to read signal for.

        returns:
            bw_idx:     Sample ID results are for.
            chromosome: Chromosome results are for.
            counts:     Array of sums for each coordinate
        """

        # Extract sum of each coordinate region
        counts = np.zeros(len(coords))

        start = 0

        for region_idx, region_coords in enumerate(coords):
            # Difference between coordinates
            region_size = region_coords[1] - region_coords[0]
            # Read signal within coordinates
            region_signal = self.signalReader(bw_idx = bw_idx,
                                              chromosome = chromosome,
                                              start_idx = region_coords[0],
                                              end_idx = region_coords[1],
                                              pad_end = False,
                                              norm_method = norm_method)

            # Count the sum for the region
            counts[region_idx] = np.sum(region_signal)
            start += region_size

        return bw_idx, chromosome, counts

    def calculateCoordCounts(self, norm_method, chromosomes = [], sample_names = []):
        """
        Calculate counts within coordinates per sample.

        params:
            norm_method:  Normalisation method of the signal to calculate sums over.
            chromosomes:  List of chromosomes to get results for.
            sample_names: List of sample names to get results for.
        """

        if len(chromosomes) == 0:
            chromosomes = self.chromosomes

        if len(sample_names) == 0:
            sample_names = self.sample_names
            sample_ids = np.arange(len(sample_names), dtype = np.uint16)
        else:
            # Get indexes matching sample names
            sample_ids = np.array([np.where(self.sample_names == s)[0][0] for s in sample_names], 
                                   dtype = np.uint16)

        counts = {}

        if self.n_cores > 1:
            with ProcessPoolExecutor(self.n_cores) as executor:
                processes = []

                for chrom in chromosomes:
                    counts[chrom] = {}
                    chrom_coords = self.coords_df.loc[self.coords_df["chrom"] == chrom][["start", "end"]].to_numpy()

                    for bw_idx in sample_ids:
                        processes.append(executor.submit(self.sampleCoordCounts,
                                                         bw_idx = bw_idx,
                                                         chromosome = chrom,
                                                         coords = chrom_coords,
                                                         norm_method = norm_method))

                for process in as_completed(processes):
                    bw_idx, chrom, results = process.result()
                    counts[chrom][bw_idx] = results

                if self.checkParallelErrors(processes):
                    raise RuntimeError("calculateCoordCounts failed to complete. To debug, see trace above.")

        else:
            for chrom in chromosomes:
                counts[chrom] = {}
                chrom_coords = self.coords_df.loc[self.coords_df["chrom"] == chrom][["start", "end"]].to_numpy()

                for bw_idx in sample_ids:
                    bw_idx, chrom, results = self.sampleCoordCounts(bw_idx = bw_idx,
                                                                    chromosome = chrom,
                                                                    coords = chrom_coords,
                                                                    norm_method = norm_method)
                    counts[chrom][bw_idx] = results
            
        return counts

    def MAPlot(self, norm_method, plot_samples = [], chromosomes = [], n_cols = 4,
               point_transparency = 0.3, point_colour = "grey", title = "", plot_width = 8, 
               plot_height = 6, pdf_name = ""):
        """
        Create MA plot(s) to compare normalisation performance per sample.

        params:
            norm_method:         Name of normalisation method to create plots for.
            plot_samples:        List of sample names to include in the plot. If not set, all samples 
                                 will be plotted.
            chromosomes:         List of chromosomes to include results for in the plot. If not set, all 
                                 chromosomes will be included.
            n_cols:              Maximum number of columns for subplots.
            point_transparency:  Value between 0 to 1 to set the transparency of scatter plot points.
            point_colour:        Name of colour for scatter plot points.
            title:               Title to display at top of plot.
            plot_width:          Length of the plot.
            plot_height:         Height of the plot.
        """

        if norm_method not in self.norm_methods:
            raise ValueError(f'Unknown normalisation method "{norm_method}"')

        if len(plot_samples) > 0:
            plot_samples = np.array(plot_samples)

            # Get sample details for the normalisation method
            norm_df = self.bigwig_df.loc[self.bigwig_df["norm"] == norm_method]

            missing_mask = np.where(~np.isin(plot_samples, norm_df["sample"].to_numpy()))[0]
            missing_samples = plot_samples[missing_mask]
            n_missing = len(missing_samples) > 0

            if n_missing > 0:
                raise ValueError(f"Invalid sample{'s' if n_missing != 1 else ''} for {norm_method}: "
                                 f'{", ".join([f"{chr(34)}{s}{chr(34)}" for s in missing_samples])}')
            
        else:
            plot_samples = self.sample_names

        if len(chromosomes) == 0:
            chromosomes = self.chromosomes

        if pdf_name:
            if not pdf_name.endswith(".pdf"):
                pdf_name = pdf_name + ".pdf"

            # Create directory to store plots
            os.makedirs(self.output_directories["plots"], exist_ok = True)

        # Get sum of signal within each region per sample and chromosome
        chrom_counts = self.calculateCoordCounts(norm_method = norm_method,
                                                 chromosomes = chromosomes,
                                                 sample_names = plot_samples)
        # Set numerical IDs for sample names
        name_id_map = {name: i for i, name in enumerate(self.sample_names)}
        sample_ids = np.array([name_id_map[s] for s in plot_samples])
        # Combine sums across chromosomes
        combined_counts = {}

        pos_sample_ids = []
        neg_sample_ids = []

        for sample_id in sample_ids:
            sample_counts = np.concatenate([chrom_counts[chrom][sample_id] for 
                                                         chrom in chrom_counts.keys()])
            
            # Split samples by positive and negative signal
            total_counts = np.sum(sample_counts)
            
            if total_counts > 0:
                pos_sample_ids.append(sample_id)
            elif total_counts < 0:
                neg_sample_ids.append(sample_id)
            else:
                if self.verbose > 0:
                    print(f"Warning: Skipping {self.sample_names[sample_id]} as counts across "
                          f"all regions were zero")
                continue

            # Get absolute value of counts to prevent log transformation errors
            sample_counts = np.abs(sample_counts)
            combined_counts[sample_id] = sample_counts

        # Create pairs of sample IDs to plot
        pos_id_pairs = self.createPairs(pos_sample_ids)
        neg_id_pairs = self.createPairs(neg_sample_ids)
        sample_id_pairs = np.concatenate((pos_id_pairs, neg_id_pairs)).astype(np.uint16)

        n_pairs = len(sample_id_pairs)
        n_cols = min(n_pairs, int(n_cols))
        n_rows = int(np.ceil(n_pairs / n_cols))

        # Create the figure
        fig, axes = plt.subplots(n_rows, n_cols, figsize = (plot_width, plot_height),
                                 constrained_layout = True)
        axes = np.atleast_1d(axes).flatten()

        for i, pair in zip(range(n_pairs), sample_id_pairs):
            ax = axes[i]

            signal_1 = combined_counts[pair[0]] + 1
            signal_2 = combined_counts[pair[1]] + 1

            # Compute M and A
            M = np.log2(signal_1) - np.log2(signal_2)
            A = 0.5 * (np.log2(signal_1) + np.log2(signal_2))

            ax.scatter(A, M, s = 10, alpha = point_transparency, color = point_colour, 
                       edgecolor = "none")
            ax.axhline(0, color = "red", linestyle = "--")
            ax.set_title(f"{self.sample_names[pair[0]]} vs {self.sample_names[pair[1]]}", fontsize = 9)
            ax.set_xlabel("A (Mean Log Count)")
            ax.set_ylabel("M (Log Fold Change)")

        for j in range(n_pairs, len(axes)):
            # Hide unused subplots
            axes[j].set_visible(False)

        if title == "":
            title = f"{norm_method} MA Plot"
            
            if len(chromosomes) == 1:
                title += f" ({chromosomes[0]})"

        fig.suptitle(title, fontweight = "bold")

        if pdf_name:
            plt.savefig(os.path.join(self.output_directories["plots"], pdf_name),
                        format = "pdf", bbox_inches = "tight")

        plt.show()

    def wassersteinDistance(self, sample_ids, bigwig_files, chromosome, start_idx = 0, end_idx = -1,
                            norm_method = "", exclude_all_zeros = False, min_values = {}, 
                            max_values = {}):
        """
        Calculate min-max scaled Wasserstein distance between two sample signals for a region.

        params:
            sample_ids:        List of sample IDs to calculate distance for.
            bigwig_files:      List of bigWig file paths with order matching sample IDs.
            chromosome:        Name of chromosome to read signal for.
            start_idx:         Start base-pair position (zero indexed).
            end_idx:           End base-pair position (zero indexed).
            norm_method:       Normalisation method to read signal for.
            exclude_all_zeros: Whether to ignore any pair of samples in which one or both signals 
                               were all zero.
            min_values:        Optionally pre-set the minimum value per sample pair to min-max scale by.
            max_values:        Optionally pre-set the maximum value per sample pair to min-max scale by.
        """

        # Create pairs of sample IDs
        sample_id_pairs = self.createPairs(sample_ids)

        signals = {}
        distances = {}

        for sample_id in sample_ids:
            # Read absolute value of signal within region for each sample
            sample_name = self.sample_names[sample_id]
            bw_file = bigwig_files[sample_name]
            signals[sample_id] = np.abs(self.extractChunkSignal(bigwig_file = bw_file,
                                                                chromosome = chromosome,
                                                                start_idx = start_idx,
                                                                end_idx = end_idx,
                                                                sample_name = sample_name,
                                                                pad_end = True))

        for pair in sample_id_pairs:
            sample_id_1 = int(pair[0])
            sample_id_2 = int(pair[1])

            # Get signals for the pair
            signal_1 = signals[sample_id_1]
            signal_2 = signals[sample_id_2]

            if exclude_all_zeros:
                if np.sum(signal_1) == 0:
                    # Ignore any pair in which one or both signals were all zero
                    continue
                if np.sum(signal_2) == 0:
                    continue

            # Min-max scale the signal
            if (sample_id_1, sample_id_2) in min_values:
                min_value = min_values[(sample_id_1, sample_id_2)]
            else:
                min_value = min(min(signal_1), min(signal_2))

            min_max_1 = signal_1 - min_value
            min_max_2 = signal_2 - min_value

            if (sample_id_1, sample_id_2) in max_values:
                max_value = max_values[(sample_id_1, sample_id_2)]
            else:
                max_value = max(max(min_max_1), max(min_max_2))

            if max_value > 0:
                min_max_1 /= max_value
                min_max_2 /= max_value

            # Calculate distance between the signals
            distances[(sample_id_1, sample_id_2)] = float(stats.wasserstein_distance(min_max_1, 
                                                                                     min_max_2))

        if norm_method != "":
            # Return the distance and an identifier
            return distances, norm_method
        else:
            return distances

    def addPlotStars(self, ax, x1, x2, y, p_value, bar_height = 0.001):
        """ 
        Add statistical significant stars to a plot.

        params:
            ax:         Axis to add stars to.
            x1:         x-axis start coordinate for bar.
            x2:         x-axis end coordinate for bar.
            y:          Lower y-axis coordinate for bar.
            p_value:    Significance value between 0 and 1.
            bar_height: Thickness of the bar.
        """

        barx = [x1, x1, x2, x2]
        bary = [y, y + bar_height, y + bar_height, y]
        ax.plot(barx, bary, c = "black")
        
        if p_value < 0.001:
            stars = "***"
        elif p_value < 0.01:
            stars = "**"
        elif p_value < 0.05:
            stars = "*"
        else:
            stars = "ns"

        ax.text((x1 + x2) / 2, y + bar_height + (0.1 * bar_height), stars, ha = "center", va = "bottom", 
                color = "black")

    def plotWasserstein(self, plot_samples = [], norm_methods = [], chromosomes = [], 
                        reference_norm = "", pair_merge_coords = True, use_chrom_maxs = False, 
                        exclude_all_zeros = False, log_scale = False, cmap = True, 
                        plot_type = "violin", title = "", plot_width = 6, plot_height = 4, 
                        pdf_name = ""):
        """
        Create violin or box plot of Wasserstein distance between samples over specified regions or
        coordinates.

        params:
            plot_samples:        List of sample names to include in the plot. If not set, all samples
                                 will be plotted.
            norm_methods:        List of normalisation methods to include in the plot. If not set, all 
                                 normalisation methods will be plotted.
            chromosomes:         List of chromosomes to include results for in the plot. If not set, all 
                                 chromosomes will be included.
            reference_norm:      Baseline normalisation to compare significance of results against.
            pair_merge_coords:   Set as True to combine region coordinates for sample pairs, or False to  
                                 use region coordinates combined across all samples. Not relavent if set 
                                 coordinates given.
            use_chrom_maxs:      Set as True to use chromosome-wide maximum values when mix-max scaling 
                                 signal during wasserstein distance calculation. Or set as False to use
                                 the maximum signal value for each coordinate pair. 
            exclude_all_zeros:   Whether to ignore any pair of samples in which one or both signals were
                                 all zero during wasserstein distance calculation.
            log_scale:           Whether to log transform the y-axis.
            cmap:                Can be set as a colour map for the violins / boxes, kept as True to use 
                                 the default colour map or False to use solid colours.
            plot_type:           Set as either "violin" or "box".
            title:               Title to display at top of plot.
            plot_width:          Length of the plot.
            plot_height:         Height of the plot.
            pdf_name:            To save plot to PDF, set this as a file name.
        """

        wasserstein_dists = {}

        if len(plot_samples) > 0:
            plot_samples = np.array(plot_samples)
            # Check if any sample names provided are not valid
            missing_mask = np.where(~np.isin(plot_samples, self.sample_names))[0]
            missing_samples = plot_samples[missing_mask]
            n_missing = len(missing_samples) > 0

            if n_missing > 0:
                raise ValueError(f"Invalid sample{'s' if n_missing != 1 else ''}: "
                                 f'"{", ".join(missing_samples)}"')
            
        else:
            plot_samples = self.sample_names

        n_samples = len(plot_samples)

        if n_samples < 2:
            raise ValueError(f"Not enough samples to create a Wasserstein distance plot. "
                             f"{n_samples} sample{'s' if n_samples != 1 else ''} were set, but a "
                             f"minimum of two is required to calculate pairwise distances.")

        if len(norm_methods) > 0:
            norm_methods = np.array(norm_methods)
            missing_mask = np.where(~np.isin(norm_methods, self.norm_methods))[0]
            missing_norm = norm_methods[missing_mask]
            n_missing = len(missing_norm) > 0

            if n_missing > 0:
                raise ValueError(f"Invalid normalisation method{'s' if n_missing != 1 else ''}: "
                                 f'"{", ".join(missing_norm)}"')
        else:
            norm_methods = self.norm_methods

        if reference_norm != "":
            if reference_norm not in norm_methods:
                raise ValueError(f'Unknown reference normalisation method "{reference_norm}"')

            # Add significant stars at top of plot to compare normalisation methods with the reference
            add_stars = True
        else:
            add_stars = False

        if len(chromosomes) == 0:
            chromosomes = self.chromosomes
        else:
            chromosomes = np.array(chromosomes)
            missing_mask = np.where(~np.isin(chromosomes, self.chromosomes))[0]
            missing_chroms = chromosomes[missing_mask]
            n_missing = len(missing_chroms) > 0

            if n_missing > 0:
                raise ValueError(f"Chromosome{'s were ' if n_missing != 1 else ' was'} "
                                 f"not found in the coordinates data: "
                                 f'"{", ".join(missing_chroms)}"')

        # Set numerical IDs for sample names
        name_id_map = {name: i for i, name in enumerate(self.sample_names)}
        sample_ids = np.array([name_id_map[s] for s in plot_samples])

        # Split samples into those with positive and negative signal
        pos_sample_ids = []
        neg_sample_ids = []

        for sample_id in sample_ids:
            # Find the highest value across chromosomes
            max_over_chroms = np.max(list(self.chrom_maxs[sample_id].values()))
            
            if max_over_chroms > 0:
                pos_sample_ids.append(sample_id)
            else:
                neg_sample_ids.append(sample_id)

        # Create each combination of samples with other samples of the same sign
        pos_id_pairs = self.createPairs(pos_sample_ids)
        neg_id_pairs = self.createPairs(neg_sample_ids)
        sample_id_pairs = np.concatenate((pos_id_pairs, neg_id_pairs)).astype(np.uint16)

        pair_chrom_maxs = {}

        if use_chrom_maxs:
            for chrom in chromosomes:
                pair_chrom_maxs[str(chrom)] = {}

                for pair in sample_id_pairs:
                    sample_id_1 = int(pair[0])
                    sample_id_2 = int(pair[1])

                    pair_chrom_maxs[chrom][(sample_id_1, sample_id_2)] = max(self.chrom_maxs[sample_id_1][chrom], 
                                                                             self.chrom_maxs[sample_id_2][chrom])

        if pair_merge_coords:
            if self.coords_df is None:
                pair_merge_coords = False

                if self.verbose > 1:
                    print(f"Using same coordinates per sample pair as only fixed coordinates were given")

            elif self.verbose > 1:
                print("Using coordinates merged per sample pair")

                if len(self.sample_names) == 2:
                    # Disable recalculation of merged coordinates
                    pair_merge_coords = False

        if title == "":
            title = "Min-Max Scaled Signal Distance"

            if len(chromosomes) == 1:
                title = title + f" ({chromosomes[0]})"

        if cmap is not None:
            add_colourmap = True

            try:
                cmap = plt.get_cmap(cmap)
            except:
                # Default colour map
                cmap = LinearSegmentedColormap.from_list("blue_to_red", ["#0000FF", "#FF0000"], 
                                                         N = 256)
        else:
            add_colourmap = False

        plot_type = plot_type.lower()

        if plot_type[0] == "v":
            plot_type = "violin"
        elif plot_type[0] == "b":
            plot_type = "box"
        else:
            raise ValueError(f'Invalid plot_type "{plot_type}". Options include: '
                             f'"violin" and "box"')

        if pdf_name:
            if not pdf_name.endswith(".pdf"):
                pdf_name = pdf_name + ".pdf"

            # Create directory to store plots
            os.makedirs(self.output_directories["plots"], exist_ok = True)

        if pair_merge_coords:
            if self.verbose > 0:
                print("Merging pair-specific peak coordinates")

            coords_dfs = {}

            for pair in sample_id_pairs:
                sample_id_1 = int(pair[0])
                sample_id_2 = int(pair[1])

                # Combine region coordinates for the pair
                regions_paths = self.regions_df[(self.regions_df["sample"] == self.sample_names[sample_id_1]) | 
                                                (self.regions_df["sample"] == self.sample_names[sample_id_2])]["regions"].to_numpy()
                
                coords_dfs[(sample_id_1, sample_id_2)] = self.mergeRegions(regions_paths = regions_paths,
                                                                           min_consensus = self.min_consensus,
                                                                           border_pad = self.border_pad)

        if self.verbose > 0:
            n_chroms = len(chromosomes)
            print(f"Calculating Wasserstein distances over {n_chroms} "
                  f"chromosome{'s' if n_chroms != 1 else ''}")

        max_values = {}

        if self.n_cores > 1:
            with ProcessPoolExecutor(self.n_cores) as executor:
                processes = []

                for norm_method in norm_methods:
                    wasserstein_dists[norm_method] = []
                    norm_bigwig_files = {}

                    norm_df = self.bigwig_df[self.bigwig_df["norm"] == norm_method]
                    norm_bigwig_files = dict(zip(norm_df["sample"], norm_df["bigwig"]))

                    for chrom in chromosomes:
                        if use_chrom_maxs:
                            max_values = pair_chrom_maxs[chrom]

                        if pair_merge_coords:
                            for pair in sample_id_pairs:
                                # Get chromosome coordinates combined across the two samples
                                sample_id_1 = pair[0]
                                sample_id_2 = pair[1]
                                pair_df = coords_dfs[(sample_id_1, sample_id_2)]
                                chrom_rows = pair_df[pair_df["chrom"] == chrom]
                                chrom_coords = chrom_rows[["start", "end"]].values.tolist()

                                for region_coords in chrom_coords:
                                    start = int(region_coords[0])
                                    end = int(region_coords[1])

                                    # Calculate distance for pair-specific coordinates
                                    processes.append(executor.submit(self.wassersteinDistance,
                                                                     sample_ids = pair,
                                                                     bigwig_files = norm_bigwig_files,
                                                                     chromosome = chrom,
                                                                     start_idx = start,
                                                                     end_idx = end,
                                                                     norm_method = norm_method,
                                                                     max_values = max_values,
                                                                     exclude_all_zeros = exclude_all_zeros))

                        else:
                            # Get chromosome coordinates to use for every sample pair
                            chrom_rows = self.coords_df[self.coords_df["chrom"] == chrom]
                            chrom_coords = chrom_rows[["start", "end"]].values.tolist()

                            for region_coords in chrom_coords:
                                start = int(region_coords[0])
                                end = int(region_coords[1])

                                # Calculate distances for all sample pairs using the same coordinates
                                processes.append(executor.submit(self.wassersteinDistance,
                                                                 sample_ids = sample_ids,
                                                                 bigwig_files = norm_bigwig_files,
                                                                 chromosome = chrom,
                                                                 start_idx = start,
                                                                 end_idx = end,
                                                                 norm_method = norm_method, 
                                                                 max_values = max_values,
                                                                 exclude_all_zeros = exclude_all_zeros))

                for process in as_completed(processes):
                    dist, norm_method = process.result()
                    wasserstein_dists[norm_method].extend(list(dist.values()))

                if self.checkParallelErrors(processes):
                    raise RuntimeError("plotWasserstein failed to complete. To debug, see trace above.")
        else:
            for norm_method in norm_methods:
                wasserstein_dists[norm_method] = []
                norm_bigwig_files = {}

                norm_df = self.bigwig_df[self.bigwig_df["norm"] == norm_method]
                norm_bigwig_files = dict(zip(norm_df["sample"], norm_df["bigwig"]))

                for chrom in chromosomes:
                    if use_chrom_maxs:
                        max_values = pair_chrom_maxs[chrom]

                    if pair_merge_coords:
                        for pair in sample_id_pairs:
                            sample_id_1 = pair[0]
                            sample_id_2 = pair[1]
                            pair_df = coords_dfs[(sample_id_1, sample_id_2)]
                            chrom_rows = pair_df[pair_df["chrom"] == chrom]
                            chrom_coords = chrom_rows[["start", "end"]].values.tolist()

                            for region_coords in chrom_coords:
                                start = int(region_coords[0])
                                end = int(region_coords[1])

                                # Calculate distance for pair-specific coordinates
                                dist, _ = self.wassersteinDistance(sample_ids = pair,
                                                                   bigwig_files = norm_bigwig_files,
                                                                   chromosome = chrom,
                                                                   start_idx = start,
                                                                   end_idx = end,
                                                                   norm_method = norm_method,
                                                                   max_values = max_values,
                                                                   exclude_all_zeros = exclude_all_zeros)
                                wasserstein_dists[norm_method].extend(list(dist.values()))


                    chrom_rows = self.coords_df[self.coords_df["chrom"] == chrom]
                    chrom_coords = chrom_rows[["start", "end"]].values.tolist()

                    for region_coords in chrom_coords:
                        start = int(region_coords[0])
                        end = int(region_coords[1])

                        # Calculate distances for all sample pairs using the same coordinates
                        dist, _ = self.wassersteinDistance(sample_ids = sample_ids,
                                                           bigwig_files = norm_bigwig_files,
                                                           chromosome = chrom,
                                                           start_idx = start,
                                                           end_idx = end,
                                                           norm_method = norm_method,
                                                           max_values = max_values,
                                                           exclude_all_zeros = exclude_all_zeros)
                        wasserstein_dists[norm_method].extend(list(dist.values()))

        if self.verbose > 0:
            print(f"Plotting {plot_type} plot")

        dist_df = pd.concat([pd.DataFrame({"norm": norm, "wasserstein_distance": dists}) for 
                             norm, dists in wasserstein_dists.items()])

        if len(dist_df) == 0:
            if self.verbose > 0:
                print("Cannot create plot as no Wasserstein distances found")

        if add_stars:
            # Map between normalisation methods to plot and indexes
            norm_id_map = {norm: i for i, norm in enumerate(norm_methods)}
            # Set values to adjust sigificant star placement
            star_base_height = max(dist_df["wasserstein_distance"])
            star_increment = star_base_height * 0.1
            star_base_height += len(norm_methods) * star_increment

        norm_pairs = self.createPairs(norm_methods)

        if add_colourmap:
            # Calculate mean wasserstein distance per normalisation method
            average_dists = [np.mean(wasserstein_dists[norm_method]) for norm_method in norm_methods]
            # Scale colours between [0,1] based on distance
            cmap_norm = mcolors.Normalize(vmin = min(average_dists), vmax = max(average_dists))
            
        # Create the plot
        fig, ax = plt.subplots(figsize = (plot_width, plot_height), constrained_layout = True)

        if plot_type == "violin":
            sns.violinplot(data = dist_df, x = "norm", y = "wasserstein_distance", ax = ax, 
                           log_scale = log_scale)
        else:
            if add_colourmap:
                # Map values to colours on the colour bar
                colour_palette = {norm_method: mcolors.to_hex(cmap(cmap_norm(dist))) for 
                                  norm_method, dist in zip(norm_methods, average_dists)}
            else:
                colour_palette = "colorblind"

            sns.boxplot(data = dist_df, x = "norm", y = "wasserstein_distance", ax = ax, 
                        log_scale = log_scale, palette = colour_palette)

        t_test_p_values = []

        for norm_pair in norm_pairs:
            norm_1 = norm_pair[0]
            norm_2 = norm_pair[1]

            # Perform t-test between Wasserstein distances of two normalisation methods
            t_test = stats.ttest_ind(wasserstein_dists[norm_1], wasserstein_dists[norm_2])
            p_value = float(t_test.pvalue)
            t_test_p_values.append(p_value)

            if add_stars:
                if (norm_1 == reference_norm) or (norm_2 == reference_norm):
                    # Add significance stars between the reference normalisation and another normalisation
                    self.addPlotStars(ax = ax, x1 = norm_id_map[norm_1], x2 = norm_id_map[norm_2], 
                                      y = star_base_height - (star_increment * norm_id_map[norm_1]), 
                                      p_value = p_value)
                    
        # Format t-test p-values as a DataFrame
        t_test_csv_file = "wasserstein_t_test_p_values.csv"
        t_test_df = pd.DataFrame({"norm_method_1": norm_pairs[:,0],
                                  "norm_method_2": norm_pairs[:,1],
                                  "t_test_p_value": t_test_p_values})
        t_test_df.to_csv(os.path.join(self.output_directories["output_stats"], t_test_csv_file), 
                         header = True, index = False)

        if add_colourmap:
            if plot_type == "violin":
                for patch_idx, violin in enumerate(ax.collections):
                    # Set colour of violin
                    norm_method = norm_methods[patch_idx]
                    mean_dist = average_dists[patch_idx]
                    violin_colour = cmap(cmap_norm(mean_dist))
                    violin.set_facecolor(violin_colour)

            sm = cm.ScalarMappable(cmap = cmap, norm = cmap_norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax = ax, orientation = "vertical", fraction = 0.02, pad = 0.03)
            cbar.set_label("Mean w", rotation = 270, labelpad = 15)

        plt.xlabel("Normalisation Method")
        plt.ylabel("Wasserstein Distance (w)")
        plt.title(title, fontweight = "bold")

        if pdf_name:
            plt.savefig(os.path.join(self.output_directories["plots"], pdf_name),
                        format = "pdf", bbox_inches = "tight")

        plt.show()