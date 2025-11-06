import numpy as np
import numba
import pandas as pd
import pyBigWig
import os
import glob
import re
import gc
import gzip
import matplotlib.pyplot as plt
import subprocess
import sympy as sp
from scipy.special import expi
from scipy import stats
from concurrent.futures import ProcessPoolExecutor
from .chrom_analysis import ChromAnalysisExtended

#####################################
# Main class for bigWig normalisation
#####################################

class ZoneNorm(ChromAnalysisExtended):
    # deepTools bamCoverage normalisation
    bam_norm_methods = ["RPKM", "CPM", "BPM", "RPGC"]
    # Inbuilt normalisation methods and other transformations
    bigwig_norm_methods = ["ZEN", "Power", "Log"]

    # Allow only supported normalisation methods
    supported_norm_methods = ["No normalisation"]
    supported_norm_methods.extend(bam_norm_methods)
    supported_norm_methods.extend(bigwig_norm_methods)

    # Allow normalisation against these options
    supported_norm_stats_types = ["signal",
                                  "signal_padded_merged_zones",
                                  "signal_unpadded_merged_zones",
                                  "signal_padded_sample_zones",
                                  "signal_unpadded_sample_zones"]

    # Supported options for test distribution parameters
    param_types = ["scipy_fit", "mean_fit", "median_fit"]

    # Map between distribution names and scipy classes
    scipy_distributions = {"norm": stats.norm,
                           "laplace": stats.laplace,
                           "logistic": stats.logistic,
                           "cauchy": stats.cauchy,
                           "gumbel_l": stats.gumbel_l,
                           "gumbel_r": stats.gumbel_r}

    # By default, only test subset of recommended distributions
    default_test_distributions = ["laplace"]

    # Use these column names in region_coords for consistency
    region_coords_cols = ["chrom", "start", "end"]

    __slots__ = ("chrom_sizes_file", "interleave_sizes", "norm_method", "genome_size", "extend_reads", 
                 "filter_strand", "exclude_zero", "zone_remove_percent", "norm_stats_type", 
                 "norm_power", "deletion_size", "downsample_size", "downsample_seed", "kernel", 
                 "test_distributions", "log_transform", "zone_distribution", "zone_param_type", 
                 "best_zone_distribution", "zone_probability", "bin_size", "extend_n_bins", 
                 "merge_depth", "min_region_bps", "quality_filter", "min_different_bps")
    def __init__(self, bam_paths = [], bigwig_paths = [], chromosomes = [], chrom_sizes_file = None, 
                 interleave_sizes = True, sample_names = [], blacklist = "", norm_method = "ZEN", 
                 genome_size = None, extend_reads = False, filter_strand = False, exclude_zero = True, 
                 zone_remove_percent = 10, norm_stats_type = "signal_padded_sample_zones",
                 norm_power = None, deletion_size = 500, downsample_size = 300, downsample_seed = 0,
                 kernel = [], test_distributions = ["laplace"], log_transform = True, 
                 zone_distribution = "laplace", zone_param_type = "median_fit", 
                 zone_probability = 0.995, bin_size = 1000, extend_n_bins = 1, merge_depth = 5, 
                 min_region_bps = 35, quality_filter = True, min_different_bps = 5, 
                 n_cores = 1, analysis_name = "Analysis", verbose = 1):
        """
        Class for predicting zones of signal and normalise bigWigs or BAMs.

        params:
            bam_paths:           List of paths of BAM files of interest. This takes priority over 
                                 bigwig_paths.
            bigwig_paths:        List of paths of bigWig and/or wig files of interest.
            chromosomes:         List of chromosomes to run analysis on.
            chrom_sizes_file:    Path to file of tab separated chromosomes and sizes for creating 
                                 bigBed or converting wig to bigWig.
            interleave_sizes:    If True and using multiple cores, then process larger chromosomes 
                                 alongside smaller ones to reduce memory usage. Otherwise process 
                                 chromosomes from largest to smallest.
            sample_names:        Optionally set as a list of custom names corresponding to each file.
                                 e.g. 'cntrl_s1.bw' and 'flt3_inhibitor_s1.bw' could be set as 
                                 ["Control Sample", "Treated Sample"].
                                 This will be converted to a dictionary mapping original file names 
                                 to the provided custom names.
                                 e.g. accessing sample_names would return 
                                 {"cntrl_s1": "Control Sample", "flt3_inhibitor_s1": "Treated Sample"}.
            blacklist:           File path to blacklist file with chromosome coordinates to exclude.
            norm_method:         Signal Normalisation method to apply, i.e. "ZEN", "Power", "Log", 
                                 "No normalisation" (bigWig/BAM) or "RPKM", "CPM", "BPM", "RPGC" (BAM).
            genome_size:         Required if using "RPGC" BAM normalisation. Can be set as either an 
                                 integer, or "hg19", "grch37", "hg38", "grch38", "t2t", "chm13cat_v2",
                                 "mm9", "grcm37", "mm10", "grcm38", "mm39", "grcm39", "dm3", "dm6",
                                 "danrer10":, "grcz10", "danrer11", "grcz11", "wbcel235" or "tair10".
            extend_reads:        Whether to enable read extension when creating bigWigs from BAMs with 
                                 deepTools.
            filter_strand:       Specifies whether to separated reads by strand when creating bigWigs 
                                 from BAMs with deepTools. Setting as True will create a forward and 
                                 reverse strand bigWig per BAM (e.g. for nascent transcription). 
                                 Leaving as False will disable strand filtering. Setting this as 
                                 "forward" or "reverse" will extract the specified strand per BAM.
            exclude_zero:        Whether to ignore zeros when calculating statistics within signal 
                                 zones.
            zone_remove_percent: Set as a value greater than zero to ignore an upper and lower 
                                 percentile when calculating statistics of signal in zones.
                                 e.g. zone_remove_percent = 100 for a signal of 1,000,000 and 
                                 norm_method = "ZEN" allows the 100th upper and lower percentile to be 
                                 removed when calculating mean and standard deviation. This prevents 
                                 extreme values, such as those caused by technical artifacts in 
                                 ATAC/ChIP-seq, from biasing the normalisation.
            norm_stats_type:     Type of statistics to use for ZEN normalisation. By default, this is 
                                 "signal_padded_sample_zones", which will calculate averages and 
                                 deviations from signal within sample-specific padded zones. Other 
                                 options include: "signal_unpadded_sample_zones", which instead uses 
                                 unpadded zones, "signal_padded_merged_zones" or 
                                 "signal_unpadded_merged_zones", which use zone coordinates merged 
                                 across all samples, or "signal", which considers all signal for 
                                 statistics calculations.
            norm_power:          If norm_method is set as "Power", set this as a number to raise 
                                 signal to the power of.
            deletion_size:       Threshold number of consecutive zeros in the signal for a region to be 
                                 considered a potential deletion. 
            downsample_size:     Number of positions of signal to select when fitting distributions.
            downsample_seed:     Integer seed for random sampling of signal when fitting distributions.
            kernel:              Smoothing kernel applied to predict signal scores. Custom kernel can 
                                 be given as an array.
            test_distributions:  List of distributions to test for signal zone prediction. 
                                 Supported distributions include: "norm", "laplace", "logistic", 
                                 "cauchy", "gumbel_l", "gumbel_r".
            log_transform:       Whether to apply logarithm before transforming the smoothed signal 
                                 during distribution testing.
            zone_distribution:   The distribution to fit for signal zone prediction. Can be set as 
                                 either a distribution name, or as "infer" to automatically set the 
                                 best tested distribution.
            zone_param_type:     The distribution parameter type to fit for signal zone distribution. 
                                 Can be set as either "mean_fit", "median_fit", "scipy_fit" or "infer" 
                                 to automatically set the best tested parameter type.
            zone_probability:    Probabilty threshold (between 0 to 1) from which the signal cut off 
                                 is dervied from the distribution fitted for signal zone prediction.
            bin_size:            Number of base pairs to use in a full bin, e.g. 1,000.
            extend_n_bins:       An integer that alongside bin size, determines the amount of padding 
                                 to add to predicted zones. E.g. if bin_size = 1,000 and 
                                 extend_n_bins = 2, an unpadded zone coordinate of [5,678, 6,789] is 
                                 first rounded to [5,000, 7,000], then extended either side by 
                                 (2 * 1,000) as [3,000, 9,000].
            merge_depth:         Minimum distance between two zones for them to be combined into one.
            min_region_bps:      Minimum size of a region initially predicted to have signal, i.e. 
                                 the number of base pairs that must consecutively exceed the zone 
                                 threshold. The zone is likely to be larger than this if rounding 
                                 region coordinates to the nearest bins and adding one or more bins 
                                 worth of padding.
            quality_filter:      Set a True to filter out regions of signal within zones with 
                                 insufficient signal.
            min_different_bps:   If using the quality filter, set the minimum number of different 
                                 base pairs that need to be found consecutively for a signal to be 
                                 considered good quality.
            n_cores:             The number of cores / CPUs to use if using multiprocessing.
            analysis_name:       Custom name of folder to save results to. By default this will be 
                                 set to "Analysis".
            verbose:             Set as an integer greater than 0 to display progress messages.
        """

        # Initialise the parent class
        super().__init__(bam_paths = bam_paths,
                         bigwig_paths = bigwig_paths,
                         allow_bams = True, # Allow bigWig, wig or BAM files to be input
                         sample_names = sample_names,
                         blacklist = blacklist,
                         n_cores = n_cores,
                         analysis_name = analysis_name,
                         verbose = verbose)

        # Set output folders
        temp_dir = self.output_directories["temp"]
        results_dir = self.output_directories["results"]

        self.output_directories = {"temp": temp_dir,
                                   "results": results_dir,
                                   "signal_zones": os.path.join(temp_dir, "Signal_Zones"),
                                   "norm_signal": os.path.join(temp_dir, "Normalised_Signals"),
                                   "smooth_signal": os.path.join(temp_dir, "Smooth_Signal"),
                                   "missing_signal": os.path.join(temp_dir, "Missing_Signal"),
                                   "signal_stats": os.path.join(temp_dir, "Signal_Stats"),
                                   "dist_stats": os.path.join(temp_dir, "Distribution_Stats"),
                                   "binned_stats": os.path.join(temp_dir, "Binned_Stats"),
                                   "bigwig": os.path.join(results_dir, "BigWigs"),
                                   "bed": os.path.join(results_dir, "BED"),
                                   "bigbed": os.path.join(results_dir, "BigBED"),
                                   "output_stats": os.path.join(results_dir, "Stats"),
                                   "smooth_bigwig": os.path.join(results_dir, "Smoothed_BigWigs"),
                                   "zones_csv": os.path.join(results_dir, "Signal_Zones"),
                                   "missing_csv": os.path.join(results_dir, "Missing_Signal"),
                                   "plots": os.path.join(self.output_directories["results"], "Plots")}

        self.setInterleaveSizes(interleave_sizes)
        self.setGenomeSize(genome_size)
        self.setNormMethod(norm_method, norm_power)
        self.setExtendReads(extend_reads)
        self.setFilterStrand(filter_strand)
        self.setExcludeZero(exclude_zero)
        self.setZoneRemovePercent(zone_remove_percent)
        self.setNormStatsType(norm_stats_type)
        self.setDeletionSize(deletion_size)
        self.setKernel(kernel)

        # Set distributions based on parameter
        self.setTestDistributions(test_distributions)
        self.setLogTransform(log_transform)
        self.setDownsampleSize(downsample_size)
        self.setDownsampleSeed(downsample_seed)
        self.setZoneDistribution(zone_distribution)
        self.setZoneParamType(zone_param_type)

        # Best distribution and parameter type to be inferred
        self.best_zone_distribution = None

        self.setZoneProbability(zone_probability)
        self.setBinSize(bin_size)
        self.setExtendNBins(extend_n_bins)
        self.setMergeDepth(merge_depth)
        self.setMinDifferentBPs(min_different_bps)
        self.setMinRegionBPs(min_region_bps)
        self.setQualityFilter(quality_filter)

        # Check whether any BAM files were given that need converting to bigWig
        if len(self.bam_paths) > 0:
            self.bamToBigWig(extend_reads = self.extend_reads, filter_strand = self.filter_strand)
            
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

        kernel_size = len(self.kernel)
        if np.all(self.kernel == self.createTriangleKernel(size = kernel_size)):
            kernel_msg = f"Triangle kernel of size {kernel_size}"
        else:
            kernel_msg = f"Custom kernel of size {kernel_size}"

        if self.blacklist is None:
            blacklist_msg = "Not set"
        else:
            blacklist_msg = self.blacklist

        dist_msg = f"{self.zone_distribution.title()} with {self.zone_param_type.replace('_', ' ')}, "
        dist_msg += f"parameter type"

        if self.verbose >= 0:
            verbose_msg = "(silent)"
        elif self.verbose == 1:
            verbose_msg = "(active)"
        else:
            verbose_msg = "(debugging mode)"

        message = (f'{self.__class__.__name__} object for "{self.analysis_name}"\n'
                   f"   * Output directory: {self.getOutputDirectory()}\n"
                   f"   * {chrom_msg}\n"
                   f"   * Number of samples: {len(self.sample_ids)}\n"
                   f"   * Sample names: {', '.join(self.sample_names.values())}\n"
                   f"   * Blacklist: {blacklist_msg}\n"
                   f"   * Kernel: {kernel_msg}\n"
                   f"   * Test distributions: {', '.join(self.test_distributions).title().replace('_', ' ')}\n"
                   f"   * Set distribution: {dist_msg}\n"
                   f"   * Zone probabilty: {self.zone_probability}\n"
                   f"   * Resources: {self.n_cores} cores\n"
                   f"   * Verbose: {self.verbose} {verbose_msg}")

        return message

    def setInterleaveSizes(self, interleave_sizes):
        try:
            self.interleave_sizes = bool(interleave_sizes)
        except:
            raise ValueError("interleave_sizes must be set as either True or False")

    def setGenomeSize(self, genome_size):
        if not (genome_size is None):
            if isinstance(genome_size, str):
                genome_size = genome_size.lower()

                genome_size_map = {"hg19": 2864785220,
                                   "grch37": 2864785220,
                                   "hg38": 2913022398,
                                   "grch38": 2913022398,
                                   "t2t": 3117292070,
                                   "chm13cat_v2": 3117292070,
                                   "mm9": 2620345972,
                                   "grcm37": 2620345972,
                                   "mm10": 2652783500,
                                   "grcm38": 2652783500,
                                   "mm39": 2654621783,
                                   "grcm39": 2654621783,
                                   "dm3": 162367812,
                                   "dm6": 142573017,
                                   "danrer10": 1369631918,
                                   "grcz10": 1369631918,
                                   "danrer11": 1368780147,
                                   "grcz11": 1368780147,
                                   "wbcel235": 100286401,
                                   "tair10": 119482012}

                if genome_size in genome_size_map:
                    self.genome_size = genome_size_map[genome_size]
            else:
                try:
                    self.genome_size = max(0, int(genome_size))
                except:
                    raise ValueError("genome_size must be set as a either a supported genome name "
                                     "or specified directly as a positive integer."
                                     "\nSupported organisms include: human (hg19/hg38), mouse "
                                     "(mm9/mm10/mm39), drosophila (dm3/dm6), "
                                     "zebrafish (danrer10/danrer11), c. elegans (wbcel235) "
                                     "and a. thaliana (tair10)")
        else:
            self.genome_size = None

    def setNormMethod(self, norm_method, norm_power = None):
        if (norm_method is None) or (norm_method == "") or (norm_method.lower()[:2] == "no"):
            self.norm_method = "No normalisation"
        else:
            # Remove special characters and convert to uppercase
            self.norm_method = "".join(filter(str.isalnum, norm_method)).upper()

            if self.norm_method == "RPGC":
                if self.genome_size is None:
                    raise ValueError("genome_size must be set to use RPGC normalisation with "
                                     "deepTools bamCoverage")

            if self.norm_method == "SQUARE" or self.norm_method == "SQUARED":
                self.norm_method = "Power"
                self.norm_power = 2

            elif self.norm_method == "CUBE" or self.norm_method == "CUBED":
                self.norm_method = "Power"
                self.norm_power = 3

            elif self.norm_method == "POWER":
                try:
                    self.norm_power = float(norm_power)
                except:
                    raise ValueError("To transform signal by a power, norm_power must be specified. "
                                     "For example, to square signal, set norm_power = 2.")
                
                self.norm_method = "Power"

            elif self.norm_method == "LOG" or self.norm_method == "LOGARITHM":
                self.norm_method = "Log"

            if len(self.bam_paths) == 0:
                if self.norm_method in self.bam_norm_methods:
                    raise ValueError(f"{self.norm_method} is only supported for BAM files")

        if self.norm_method not in self.supported_norm_methods:
            raise ValueError(f'Normalisation method "{norm_method}" not supported.\n'
                             f"Set norm_method as one of the following: "
                             f'{", ".join(self.supported_norm_methods)}'
                             f", or leave as None to disable this step.")

    def setExtendReads(self, extend_reads):
        extend_reads = str(extend_reads).strip().lower()

        if extend_reads in ["false", "none", "no", "0"]:
            self.extend_reads = False
            return
        if extend_reads in ["true", "yes", "1"]:
            self.extend_reads = True
            return
        
        raise ValueError('extend_reads must be either "True" or "False", or "None" if not '
                         "being used")

    def setFilterStrand(self, filter_strand):
        if filter_strand is None:
            self.filter_strand = False
            return
        
        elif isinstance(filter_strand, str):
            filter_strand = filter_strand.strip().lower()
            if filter_strand in ["forward", "reverse"]:
                self.filter_strand = filter_strand
                return
            if filter_strand in ["false", "none", "no", "0"]:
                self.filter_strand = False
                return
            if filter_strand in ["true", "yes", "1"]:
                self.filter_strand = True
                return
            
        if isinstance(filter_strand, bool):
            self.filter_strand = filter_strand
            return
        
        raise ValueError('filter_strand must be set as either: "forward", "reverse", "True" or "False"')

    def setExcludeZero(self, exclude_zero):
        try:
            self.exclude_zero = bool(exclude_zero)
        except:
            raise ValueError("exclude_zero must be set as either True or False "
                             "for zone statistics calculations")

    def setZoneRemovePercent(self, zone_remove_percent):
        try:
            self.zone_remove_percent = max(int(zone_remove_percent), 0)
        except:
            raise ValueError("zone_remove_percent must be set as an integer value, e.g. 10")

    def setNormStatsType(self, norm_stats_type):
        if norm_stats_type not in self.supported_norm_stats_types:
            raise ValueError(f"Invalid norm_stats_type {norm_stats_type}."
                             f"norm_stats_type must be set as either "
                             f'{", ".join(self.supported_norm_stats_types)}')
        else:
            self.norm_stats_type = norm_stats_type

    def setDeletionSize(self, deletion_size):
        try:
            # Set limit between 10 and 10000 base pairs
            self.deletion_size = min(10000, max(int(deletion_size), 10))
        except:
            raise ValueError("deletion_size must be set as an integer between 10 and 10000")

    def setKernel(self, kernel = []):
        try:
            # Allow no kernel to be used
            if kernel is not None:
                if len(kernel) == 0:
                    # Set default smoothing kernel
                    self.kernel = self.createTriangleKernel(size = 301)
                else:
                    # Use custom kernel
                    self.kernel = np.array(kernel, dtype = np.float32)
        except:
            raise ValueError("kernel must be set as a list/array."
                             "The default kernel will be used if this parameter is set as "
                             "an empty list, i.e. [].")


    def setTestDistributions(self, test_distributions):
        """ 
        Sets or updates the distributions to test for zone prediction.
        
        params:
            test_distributions: List of new distributions to test.
        """

        # Convert to list if not already
        test_distributions = list(test_distributions)
        supported_distributions = self.getSupportedDistributions()

        if len(test_distributions) > 0:
            # Check if given distributions are supported
            if np.all(np.isin(np.array(test_distributions), np.array(supported_distributions))):
                self.test_distributions = test_distributions

            else:
                raise ValueError(f"test_distributions must be set as a list of distribution names, "
                                 f"or a dictionary of SciPy distributions.\n"
                                 f"Supported distributions include: "
                                 f'{", ".join(supported_distributions)}')
        else:
            self.test_distributions = self.default_test_distributions

    def setLogTransform(self, log_transform, verbose = 0):
        try:
            self.log_transform = bool(log_transform)
        except:
            raise ValueError("log_transform must be set as True or False.")
            
        if verbose > 0:
            print(f"Updating log tranformation as {log_transform}")

    def setDownsampleSize(self, downsample_size, verbose = 0):
        try:
            # Set downsample size
            downsample_size = int(downsample_size)

            if downsample_size < 100:
                print("downsample_size is too small. Setting as the default of 300.")
                downsample_size = 300
            
            self.downsample_size = downsample_size
        except:
            raise ValueError("downsample_size must be set as an integer 100 or more")

        if verbose > 0:
            print(f"Updating downsample size as {downsample_size}")

    def setDownsampleSeed(self, downsample_seed):
        try:
            # Set downsample seed
            self.downsample_seed = max(0, int(downsample_seed))
        except:
            raise ValueError("downsample_seed must be set as an integer")

    def setZoneDistribution(self, zone_distribution):
        supported_dists = self.getSupportedDistributions()

        if zone_distribution in supported_dists:
            self.zone_distribution = zone_distribution
        else:
            zone_distribution = zone_distribution.lower()
            if zone_distribution in ["best", "infer"]:
                self.zone_distribution = "infer"
            else:
                raise ValueError(f'zone_distribution must be set as either "infer" '
                                 f"or the name of a distribution in test_distributions.\n"
                                 f"Current options include: ",
                                 f"{', '.join(supported_dists)}")

        if zone_distribution not in self.test_distributions:
            if self.verbose > 0:
                print(f"Adding {zone_distribution} to test distributions")

            self.test_distributions.append(zone_distribution)

    def setZoneParamType(self, zone_param_type):
        zone_param_type = zone_param_type.lower().removesuffix("_fit")

        if zone_param_type in [p.removesuffix("_fit") for p in self.param_types]:
            self.zone_param_type = zone_param_type + "_fit"
        else:
            zone_param_type = zone_param_type.lower()
            if zone_param_type in ["best", "infer"]:
                self.zone_param_type = "infer"
            else:
                raise ValueError(f'zone_param_type must be set as either "infer", '
                                 f'"scipy_fit", "mean_fit", or "median_fit"')

    def setZoneProbability(self, zone_probability):
        try:
            self.zone_probability = min(1, max(zone_probability, 0))
        except:
            raise ValueError("zone_probability must be set as a probability between 0 to 1")

    def setBinSize(self, bin_size = None):
        try:
            # Set limit between 100 and 10000 base pairs
            self.bin_size = min(10000, max(int(bin_size), 100))
        except:
            raise ValueError("bin_size must be set as an integer between 100 and 10000")
    
    def setExtendNBins(self, extend_n_bins):
        try:
            self.extend_n_bins = max(1, min(int(extend_n_bins), 10))
        except:
            raise ValueError("extend_n_bins must be set as an integer between 1 to 10")

    def setMergeDepth(self, merge_depth):
        try:
            self.merge_depth = max(0, min(int(merge_depth), 10000))
        except:
            raise ValueError("merge_depth must be set as an integer between 0 to 10000")

    def setMinDifferentBPs(self, min_different_bps):
        try:
            self.min_different_bps = max(1, min(int(min_different_bps), 100))
        except:
            raise ValueError("min_different_bps must be set as an integer between 1 to 1000")

    def setMinRegionBPs(self, min_region_bps):
        try:
            self.min_region_bps = max(10, min(int(min_region_bps), 10000))
        except:
            raise ValueError("min_region_bps must be set as an integer between 10 and 10000")

    def setQualityFilter(self, quality_filter):
        if isinstance(quality_filter, bool):
            self.quality_filter = quality_filter
        else:
            quality_filter = str(quality_filter).lower()
            if (quality_filter[0] == "f") or (quality_filter[:2] == "no"):
                self.quality_filter = False
            else:
                # Default to applying a quality filter
                self.quality_filter = True


    def runBamCoverage(self, bam_file, bigwig_file, subprocess_cores, extend_reads, norm_method, 
                       scale_factor = 1, effective_genome_size = 0, filter_strand = "", 
                       bin_size = 1, sample_name = None, subprocess_verbose = None):
        """
        Create a bigWig from a BAM using deepTools bamCoverage.
        
        params:
            bam_file:              BAM file to process.
            bigwig_file:           Name of output file.
            subprocess_cores:      Number of cores to use for multiprocessing.
            extend_reads:          Set as False to disable read extension.
            norm_method:           Set as either "No Normalisation" or a supported deepTools 
                                   normalisation method.
            scale_factor:          Value to scale reads by.
            effective_genome_size: Size of mapped genome (required if using RPGC normalisation).
            filter_strand:         Leave empty to disable strand filtering, or set as "forward" 
                                   or "reverse". 
                                   Filtering by reverse signal will create negative signal (if 
                                   scale_factor is not set).
            bin_size:              Base pair resolution of coverage track.
            sample_name:           Sample name can be provided for printing progress message.
            subprocess_verbose:    Verbose of stdout when running subprocess.
        """

        # Set the bamCoverage command
        command = f"bamCoverage -b {bam_file} -o {bigwig_file} -bs {bin_size} "
        command += f"--numberOfProcessors {subprocess_cores}"

        # Add additional parameters
        if extend_reads:
            command += " --extendReads"

        if norm_method != "No Normalisation":
            command += f" --normalizeUsing {norm_method}"

        if effective_genome_size > 0:
            command += f" --effectiveGenomeSize {effective_genome_size}"

        if filter_strand:
            command += f" --filterRNAstrand {filter_strand}"

            if filter_strand == "reverse":
                # Create negative signal for reverse stand
                if scale_factor == 1:
                    scale_factor = -1

        if scale_factor != 1:
            command += f" --scaleFactor {scale_factor}"

        if self.verbose > 0:
            message = f"Creating {filter_strand}{' strand ' if filter_strand else ''}"
            message += f"bigWig with deepTools"

            if sample_name is not None:
                message += f" for {sample_name}"
            
            print(message)

        # Run bamCoverage
        subprocess.run(command, shell = True, stdout = subprocess_verbose)


    def bamToBigWig(self, bam_paths = [], bin_size = 1, norm_method = "", subprocess_cores = 6,
                    effective_genome_size = None, extend_reads = None, filter_strand = False, 
                    replace_existing = False):
        """
        Run deepTools to create bigWigs from one or more BAM files.

        params:
            bam_paths:             List of BAM file paths to create bigWigs for
            bin_size:              Base pair resolution of coverage track.
            norm_method:           Set as either "No Normalisation" or a supported deepTools 
                                   normalisation method.
            subprocess_cores:      Number of cores to use per bigWig.
            effective_genome_size: Size of mapped genome (required if using RPGC normalisation).
            extend_reads:          Set as False to disable read extension.
            filter_strand:         Set as False to disable strand filtering, True to use both 
                                   strands, or set as "forward" or "reverse".
            replace_existing:      Whether to overwrite previously created files.
        """

        if len(norm_method) == 0:
            norm_method = self.norm_method

            if norm_method not in self.bam_norm_methods:
                # Normalisation will be applied later
                norm_method = "No Normalisation"

        elif norm_method not in self.bam_norm_methods:
            raise ValueError(f"{self.norm_method} is not supported for deepTools bamCoverage")

        if norm_method != "No Normalisation":
            output_prefix = f"_{norm_method}"
            output_directory = os.path.join(self.output_directories["bigwig"], norm_method)
        else:
            output_prefix = ""
            output_directory = os.path.join(self.output_directories["bigwig"], "No_Normalisation")

        os.makedirs(output_directory, exist_ok = True)

        if effective_genome_size is None:
            if not (self.genome_size is None):
                effective_genome_size = self.genome_size
            else:
                effective_genome_size = 0

        if filter_strand:
            filter_strand_dict = {"forward": "Pos",
                                  "reverse": "Neg"}
            if filter_strand in filter_strand_dict.keys():
                # Subset dictionary for specific strand
                filter_strand_dict = {filter_strand: filter_strand_dict[filter_strand]}
            elif filter_strand == True:
                pass
            else:
                raise ValueError(f'Unknown value "{filter_strand}" for filter_strand')
        else:
            filter_strand_dict = {}
            strand = ""

        if self.verbose > 0:
            subprocess_verbose = None
        else:
            # Disable any printing
            subprocess_verbose = subprocess.DEVNULL

        if len(bam_paths) == 0:
            set_attributes = True
            bam_paths = self.bam_paths
            bigwig_paths = []
            update_sample_names = []
        else:
            # Cannot update attributes if custom list was given as it may not map to samples
            set_attributes = False

        # Cap cores to not exceed total avaliable
        subprocess_cores = min(int(subprocess_cores), self.n_cores)

        # Record file names of BAMs to process and bigWigs to create
        process_bam_files = {}
        process_bigwig_files = {}
        create_n_bigwigs = 0

        for bam_file in bam_paths:
            # Use BAM file path to set bigWig file names
            sample_name = bam_file.split(os.sep)[-1].replace(".bam", "")

            if len(filter_strand_dict) > 0:
                bigwig_files = []
                for strand in filter_strand_dict:
                    # Set to create strand specific bigWigs
                    bigwig_files.append(os.path.join(output_directory, 
                                                     f"{sample_name}_{filter_strand_dict[strand]}{output_prefix}.bw"))
            else:
                # Single bigWig
                bigwig_files = [os.path.join(output_directory, f"{sample_name}{output_prefix}.bw")]
            
            if set_attributes:
                bigwig_paths.extend(bigwig_files)
                update_sample_names.append(sample_name)

            if (not replace_existing) and (np.all([os.path.exists(f) for f in bigwig_files])):
                try:
                    for f in bigwig_files:
                        # Check if the bigWig can be opened without error
                        with pyBigWig.open(f) as bw:
                            pass
                    # Skip creating a bigWig that already exists
                    continue

                except (OSError, RuntimeError):
                    # Allow the bigWig(s) to be recreated if an error occured
                    pass

            # Record files to process
            process_bam_files[sample_name] = bam_file
            process_bigwig_files[sample_name] = bigwig_files
            create_n_bigwigs += len(bigwig_files)

        if create_n_bigwigs:
            if self.n_cores > 1:
                # Check if underutilising cores
                if create_n_bigwigs * subprocess_cores < self.n_cores:
                    # Increate cores per process
                    subprocess_cores = self.n_cores // create_n_bigwigs

                executor = ProcessPoolExecutor(self.n_cores // subprocess_cores)
                processes = []

            for sample_name, bigwig_files in process_bigwig_files.items():
                bam_file = process_bam_files[sample_name]

                for i, bigwig_file in enumerate(bigwig_files):
                    if len(filter_strand_dict) > 0:
                        # Get forward or reverse strand
                        strand = list(filter_strand_dict.keys())[i]

                    if self.n_cores > 1:
                        processes.append(executor.submit(self.runBamCoverage,
                                                         bam_file = bam_file,
                                                         bigwig_file = bigwig_file,
                                                         bin_size = bin_size,
                                                         subprocess_cores = subprocess_cores,
                                                         extend_reads = extend_reads,
                                                         filter_strand = strand,
                                                         norm_method = norm_method,
                                                         effective_genome_size = effective_genome_size,
                                                         sample_name = sample_name,
                                                         subprocess_verbose = subprocess_verbose))
                        
                    else:
                        self.runBamCoverage(bam_file = bam_file,
                                            bigwig_file = bigwig_file,
                                            bin_size = bin_size,
                                            subprocess_cores = subprocess_cores,
                                            extend_reads = extend_reads,
                                            filter_strand = strand,
                                            norm_method = norm_method,
                                            effective_genome_size = effective_genome_size,
                                            sample_name = sample_name,
                                            subprocess_verbose = subprocess_verbose)

            if self.n_cores > 1:
                if self.checkParallelErrors(processes):
                    return None

        if set_attributes:
            self.bigwig_paths = np.array(bigwig_paths)

            if filter_strand:
                # Update sample names if split into positive and negative strands
                for sample_name in update_sample_names:
                    if sample_name in self.sample_names:
                        # Remove original name from sample names
                        custom_sample_name = self.sample_names[sample_name]
                        del self.sample_names[sample_name]
                    else:
                        custom_sample_name = sample_name

                    # Add positive and negative varients
                    self.sample_names[f"{sample_name}_Pos"] = f"{custom_sample_name}_Pos"
                    self.sample_names[f"{sample_name}_Neg"] = f"{custom_sample_name}_Neg"

                self.file_sample_names = [s.replace("_", "-").replace(".", "-") for 
                                          s in self.sample_names.keys()]
                # Update IDs
                self.sample_ids = np.arange(len(self.sample_names), dtype = np.uint16)

    @staticmethod
    def wigToBigWig(wig_file, chrom_sizes_file, subprocess_verbose = None):
        """
        Convert wig file to bigWig

        params:
            wig_file:           File path of wig file to convert.
            chrom_sizes_file:   Path to file of tab separated chromosomes and sizes.
            subprocess_verbose: Verbose of stdout when running subprocess.
        """

        if not os.path.exists(chrom_sizes_file):
            raise FileNotFoundError(f'Could not find chrom_sizes_file "{chrom_sizes_file}"')

        sample_name = wig_file.split(os.sep)[-1].replace(".wig", "")
        command = f"wigToBigWig {wig_file} {chrom_sizes_file} {sample_name}.bw"
        subprocess.run(command, shell = True, stdout = subprocess_verbose)

    def openDefaultSignal(self, bw_idx, chromosome, signal_type = ""):
        chrom_signal = np.array([])

        if signal_type and signal_type.lower()[0] == "s":
            signal_type = "Smoothed"
            directory = self.output_directories["smooth_signal"]
            file_prefix = os.path.join(chromosome, f"smooth-signal")
        else:
            signal_type = "Normalised"
            directory = self.output_directories["norm_signal"]
            file_prefix = os.path.join(self.norm_method, chromosome, f"{self.norm_method}-norm-signal")

        try:
            # Open numpy array of normalised signal for the chromosome
            chrom_signal = self.readArrayFile(bw_idx = bw_idx,
                                              directory = directory,
                                              file_prefix = file_prefix)[1]
        except Exception as e:
            if self.verbose > 0:
                sample_name = self.getSampleNames(return_custom = True)[bw_idx]
                print(f"Warning: {signal_type} signal not found for {chromosome} {sample_name} "
                      f"due to exception:\n{e}")
                
        return chrom_signal

    @staticmethod
    def createTriangleKernel(size = 301):
        """ 
        Create a triangular smoothing kernel so that larger weight is given to the centre of the 
        signal.
        
        params:
            size: Kernel length in base pairs.

        """

        left_size = np.ceil(size / 2)
        right_size = np.floor(size / 2)

        kernel = np.concatenate((np.arange(1, left_size + 1), 
                                 np.flip(np.arange(1, right_size + 1)))).astype(np.float32)
        kernel = kernel / np.sum(kernel)

        return kernel

    def getTestDistributions(self):
        """
        Returns distributions set to test for separating signal from background.
        """
        return self.test_distributions

    @classmethod
    def getSupportedDistributions(cls):
        """
        Return a list of distributions that can be used for signal zone prediction.
        """

        return list(cls.scipy_distributions.keys())

    def getZoneProbability(self):
        """
        Returns the set zone probability.
        """

        return self.zone_probability

    def getBlacklist(self, blacklist_file = None, chromosome = None):
        """ 
        Opens the blacklist regions and returns them as a DataFrame.
        
        params:
            blacklist_file: File path to a zipped BED file of blacklist regions to exclude, 
                            e.g. "hg38-blacklist.v2.bed.gz".
            chromosome:     Setting this as a chromosome will subset blacklist regions to those within 
                            the chromosome.
        """

        if blacklist_file is None:
            if self.blacklist is None:
                print("No blacklist set")
                return None
            else:
                blacklist_file = self.blacklist

        blacklist_df = pd.read_csv(blacklist_file, sep = "\t", header = None, 
                                   names = ["chrom", "start", "end", "class"])

        if chromosome is not None:
            # Subset regions
            blacklist_df = blacklist_df.loc[blacklist_df["chrom"] == chromosome]

        return blacklist_df

    def getStats(self, stats_type, file_name = "", lock = None, replace_existing = False):
        """
        Read statistics about the signal or distribution fitting from CSV.
        
        params:
            stats_type:       Either 'signal' or 'distribution' statistics.
            file_name:        If provided, read directly from a file.
            lock:             Can be set as a threading/multiprocessing lock to prevent 
                              parallel access to the signal statistics.
            replace_existing: Whether to overwrite previously created files.
        """

        stats_type = stats_type.lower()
        get_mean_stats = False

        if stats_type.startswith("s"):
            stats_type = "signal"
        elif stats_type.startswith("d"):
            if stats_type.endswith("mean"):
                get_mean_stats = True
            stats_type = "distribution"
        else:
            raise ValueError('stats_type must be set as either "signal" or "distribution"')

        if len(file_name) == 0:
            if stats_type == "signal":
                # Set file as combined signal statistics
                file_name = os.path.join(self.output_directories["output_stats"],
                                         "combined-signal-stats.csv")
            else:
                if get_mean_stats:
                    # Set file as combined distribution statistics
                    file_name = os.path.join(self.output_directories["output_stats"],
                                             "mean-distribution-stats.csv")
                else:
                    # Set file as combined distribution statistics
                    file_name = os.path.join(self.output_directories["output_stats"],
                                             "combined-distribution-stats.csv")

        if not replace_existing and os.path.isfile(file_name):
            if lock is not None:
                with lock:
                    # Open the statistics
                    stats_df = pd.read_csv(file_name)
            else:
                stats_df = pd.read_csv(file_name)

            if not get_mean_stats:
                # Ensure chromosome column is correct type
                stats_df["chrom"] = stats_df["chrom"].astype(str)
                
            if not stats_df.empty:
                return stats_df

        # Create the statistics
        if get_mean_stats:
            self.inferBestDistribution()
            if not os.path.exists(file_name):
                raise ValueError(f'Cannot find mean distribution statistics file "{file_name}"')

            # Get average statistics
            stats_df = pd.read_csv(file_name)
            
        else:
            # Get statistics per sample
            stats_df = self.combineStats(stats_type = stats_type, file_name = file_name, return_df = True)
            stats_df["chrom"] = stats_df["chrom"].astype(str)

        return stats_df
    
    def getSignalStats(self, file_name = "", lock = None, replace_existing = False):
        """
        Read signal statistics from CSV.

        params:
            file_name:        If provided, read directly from a file.
            lock:             Can be set as a threading/multiprocessing lock to prevent 
                              parallel access to the signal statistics.
            replace_existing: Whether to overwrite previously created files.
        """

        return self.getStats(stats_type = "signal", file_name = file_name, lock = lock, 
                             replace_existing = replace_existing)

    def getDistributionStats(self, file_name = "", get_mean_stats = False, lock = None, 
                             replace_existing = False):
        """
        Read distribution fitting statistics from CSV.

        params:
            file_name:        If provided, read directly from a file.
            lock:             Can be set as a threading/multiprocessing lock to prevent 
                              parallel access to the signal statistics.
            get_mean_stats:   
            replace_existing: Whether to overwrite previously created files.
        """

        if get_mean_stats:
            return self.getStats(stats_type = "distribution_mean", file_name = file_name, lock = lock, 
                                 replace_existing = replace_existing)
        else:
            return self.getStats(stats_type = "distribution", file_name = file_name, lock = lock, 
                                 replace_existing = replace_existing)

    def getChromFragMean(self, chromosome, sample_ids = [], signal_type = "signal_non_zero", 
                         file_name = ""):
        """
        Get the minimum and mean non-zero signal across a chromosome for a set of samples.
        
        params:
            chromosome:  The chromosome to get results for, e.g. "chr1".
            sample_ids:  The numerical IDs of samples to get results for. 
                         If left empty, then results will be returned across all samples.
            signal_type: The name of the signal to get results for. Defaults to the original 
                         signal excluding zeros.
            file_name:   File path to the signal statistics CSV. If not set, uses the default path.
        """

        if len(sample_ids) == 0:
            # Use all samples
            sample_ids = self.sample_ids

        # Convert from sample ID to name
        sample_names = np.array(self.getSampleNames(return_custom = True))[sample_ids]

        # Open the combined statistics calculated for the signal
        stats_df = self.getSignalStats(file_name = file_name)

        for column in ["sample", "chrom", "signal_type", "mean"]:
            if column not in stats_df.columns:
                raise ValueError(f'Combined signal statistics is missing column "{column}"')

        # Extract rows for the specified signal and samples
        stats_df = stats_df.loc[(stats_df["signal_type"] == signal_type) &
                                (stats_df["sample"].isin(sample_names))]

        if len(stats_df) == 0:
            # Force recreation of file if missing rows
            stats_df = self.getSignalStats(file_name = file_name, replace_existing = True)
            stats_df = stats_df.loc[(stats_df["signal_type"] == signal_type) &
                                    (stats_df["sample"].isin(sample_names))]
            
            if len(stats_df) == 0:
                raise ValueError(f"No rows found in combined signal statistics with signal "
                                 f'type "{signal_type}"')

        try:
            # Calculate average mean across chromosomes per sample
            global_mean_df = pd.DataFrame(stats_df.groupby(["sample"])["mean"].mean())
            # Reorder to ensure samples match bigWig IDs
            global_mean_df = global_mean_df.reindex(sample_names).reset_index()
            global_means = np.array(global_mean_df["mean"])

            # Extract rows for the chromosome
            stats_df = stats_df.loc[stats_df["chrom"] == chromosome]
            # Keep only sample names, estimated fragment size and non-zero mean
            stats_df = stats_df[["sample", "fragmentEstimate", "mean"]]

        except Exception as e:
            print(f"Could not extract global and {chromosome} statistics from "
                  f'"{file_name}" due to exception: ')
            print(e)
            return None

        row_sample_names = np.array(stats_df["sample"])
        chrom_fragments = np.array(stats_df["fragmentEstimate"])
        chrom_means = np.array(stats_df["mean"])

        # Check if any samples missing in the stats data
        missing_samples = np.setdiff1d(sample_names, row_sample_names)

        if len(missing_samples) > 0:
            raise ValueError(f"Cannot find {len(missing_samples)} sample(s) for {chromosome} in "
                             f'combined-signal-stats.csv: {", ".join(missing_samples)}')

        # Sort rows by sample to ensure stats align
        stats_df["sample"] = pd.Categorical(stats_df["sample"], categories = sample_names,
                                            ordered = True)
        stats_df = stats_df.sort_values("sample")

        # Extract statistics per sample
        chrom_fragments = np.array(stats_df["fragmentEstimate"])
        chrom_means = np.array(stats_df["mean"])

        return chrom_fragments, chrom_means, global_means

    def getOriginalSignal(self, sample, chromosome):
        """
        Read original bigWig signal for a specific sample and chromsome into an array.

        params:
            sample_name: Name of sample to get signal for.
            chromosome: Chromosome to get signal for.
        """

        custom_sample_names = np.array(self.getSampleNames(return_custom = True))

        if isinstance(sample, str):
            bw_idx = self.sampleToIndex(sample)
        elif isinstance(sample, (int, np.integer)):
            bw_idx = sample
            sample = custom_sample_names[bw_idx]
        else:
            raise ValueError("Invalid type for sample")

        if self.verbose > 0:
            print(f"Reading {chromosome}:0-{self.chrom_sizes[chromosome]} for {sample}")

        # Read whole chromosome signal from bigWig
        signal = self.signalReader(bw_idx = bw_idx,
                                   start_idx = 0,
                                   end_idx = -1,
                                   chromosome = chromosome,
                                   pad_end = True)
        
        return signal

    def getSmoothedSignal(self, sample, chromosome):
        """
        Read smoothed signal for a specific sample and chromsome into an array.

        params:
            sample_name: Name of sample to get signal for.
            chromosome: Chromosome to get signal for.
        """

        custom_sample_names = np.array(self.getSampleNames(return_custom = True))

        if isinstance(sample, str):
            bw_idx = self.sampleToIndex(sample)
        elif isinstance(sample, (int, np.integer)):
            bw_idx = sample
            sample = custom_sample_names[bw_idx]
        else:
            raise ValueError("Invalid type for sample")

        if self.verbose > 0:
            print(f"Reading {chromosome}:0-{self.chrom_sizes[chromosome]} for {sample}")

        # Read whole chromosome signal from array
        smooth_signal = self.readArrayFile(bw_idx = bw_idx,
                                           directory = self.output_directories["smooth_signal"],
                                           file_prefix = os.path.join(chromosome, "smooth-signal"))[1]

        return smooth_signal

    def getNormSignal(self, sample, chromosome):
        """
        Read normalised bigWig signal for a specific sample and chromsome into an array.

        params:
            sample_name: Name of sample to get signal for.
            chromosome: Chromosome to get signal for.
        """

        custom_sample_names = np.array(self.getSampleNames(return_custom = True))

        if isinstance(sample, str):
            bw_idx = self.sampleToIndex(sample)
        elif isinstance(sample, (int, np.integer)):
            bw_idx = sample
            sample = custom_sample_names[bw_idx]
        else:
            raise ValueError("Invalid type for sample")

        if self.verbose > 0:
            print(f"Reading {chromosome}:0-{self.chrom_sizes[chromosome]} for {sample}")

        # Read whole chromosome signal from bigWig
        norm_signal = self.signalReader(bw_idx = bw_idx,
                                        start_idx = 0,
                                        end_idx = -1,
                                        chromosome = chromosome,
                                        pad_end = True,
                                        read_normalised = True)
        
        return norm_signal

    def getSampleZones(self, chromosome, get_unpadded = True, get_padded = True):
        """
        Open signal zones coordinates for each sample. 
        These must first have been saved to a gzipped file.

        params:
            chromosome:   The chromosome to get signal zones for, e.g. "chr1"
            get_unpadded: Set as True to include unpadded signal zones.
            get_padded:   Set as True to include padded signal zones.
            
        returns:
            chrom_signal_zones: Dictionary of unpadded and/or padded signal zone coordinates 
                                for the sample for a specific chromosome.
        """

        custom_sample_names = np.array(self.getSampleNames(return_custom = True))

        all_zones = {}
        zone_types = []

        if get_unpadded:
            zone_types.append("unpadded")
        if get_padded:
            zone_types.append("padded")

        for zone_type in zone_types:
            sample_zones = {}

            for bw_idx in self.sample_ids:
                sample_name = str(custom_sample_names[bw_idx])
                zone_file = os.path.join(self.output_directories["signal_zones"], 
                                         chromosome, 
                                         zone_type.title(), 
                                         f"{zone_type}-zones_{self.file_sample_names[bw_idx]}.npy.gz")

                if os.path.isfile(zone_file):
                    # Record zone coordinates for the sample
                    zones = self.loadGzip(file_name = zone_file)
                    sample_zones[sample_name] = zones

            all_zones[zone_type] = sample_zones

        return all_zones
    
    def getMergedZones(self, chromosome, get_unpadded = True, get_padded = True):
        """ 
        Open signal zones merged across all samples. 
        These must first have been saved to a gzipped file.
        
        params:
            chromosome:   The chromosome to get merged zones for, e.g. "chr1"
            get_unpadded: Set as True to include unpadded signal zones.
            get_padded:   Set as True to include padded signal zones.

        returns:
            merged_zones: Dictionary of unpadded and/or padded signal zones for the chromosome.
        """

        merged_zones = {}
        zone_types = []

        if get_unpadded:
            zone_types.append("unpadded")
        if get_padded:
            zone_types.append("padded")

        for zone_type in zone_types:
            zone_file = os.path.join(self.output_directories["signal_zones"], 
                                     chromosome, zone_type.title(), f"{zone_type}_merged_zones.npy.gz")
            
            if os.path.isfile(zone_file):
                try:
                    # Open zones overlapped across samples
                    merged_zones[zone_type] = self.loadGzip(file_name = zone_file)
                except Exception as e:
                    print(e)

        return merged_zones

    def signalReader(self, bw_idx, chromosome, start_idx = 0, end_idx = -1, pad_end = False, 
                     read_normalised = False, dtype = np.float32, verbose = 0):
        """
        Read bigWig signal into a numpy array.

        params:
            bw_idx:          Index of sample to read signal for.
            chromosome:      Chromosome to read signal for.
            start_idx:       Start base-pair position (zero indexed).
            end_idx:         End base-pair position (zero indexed).
            pad_end:         Whether to add zeros to end to ensure array size is consistent 
                             across samples for a chromosome.
            read_normalised: If True, then attempt to read normalised signal.
            dtype:           Data type to store the signal as.
            verbose:         Set as a number > 0 to print progress.
        """

        sample_name = self.getSampleNames(return_custom = True)[bw_idx]

        if read_normalised:
            # Open the normalised bigWig for the sample
            bigwig_file = os.path.join(self.output_directories["bigwig"], self.norm_method,
                                       f"{sample_name}_{self.norm_method}.bw")
        else:
            # Open the original bigWig for the sample
            bigwig_file = self.bigwig_paths[bw_idx]

        # Split signal reading into separate function as these remain the same, but prior steps 
        # can be overwritten
        return self.extractChunkSignal(bw_idx = bw_idx,
                                       bigwig_file = bigwig_file,
                                       sample_name = sample_name,
                                       chromosome = chromosome,
                                       start_idx = start_idx,
                                       end_idx = end_idx,
                                       pad_end = pad_end,
                                       dtype = dtype,
                                       verbose = verbose)

    def createSignalMask(self, bw_idx, chromosome, mask_size, missing_regions = None):
        """
        Create a boolean mask to hide regions of missing signal.

        params:
            bw_idx:          Index of sample to create mask for.
            chromosome:      Chromosome to create mask for.
            mask_size:       Signal length.
            missing_regions: Coordinates of signal regions with all-zero signal.
        """

        if missing_regions is None:
            try:
                missing_regions = self.readArrayFile(bw_idx = bw_idx,
                                                     directory = self.output_directories["missing_signal"],
                                                     file_prefix = os.path.join(chromosome, 
                                                                                "missing-regions"))[1]
            except:
                if self.verbose > 0:
                    # Missing regions could be absent if the chromosome was sparse
                    sample_name = self.getSampleNames(return_custom = True)[bw_idx]
                    print(f"Warning: Could not find missing regions for {chromosome} {sample_name}")
                return None

        # Create mask of boolean values to select signal and exclude regions of zeros
        signal_mask = np.full(mask_size, True)

        for region in missing_regions:
            signal_mask[region[0]:region[1]] = False

        return signal_mask

    @staticmethod
    @numba.njit
    def numbaMaskMissing(signal, deletion_size):
        """
        Create signal mask to hide regions with gaps of zeros and return coordinates of 
        non-masked regions.
        
        params:
            signal:        Signal to mask gaps for.
            deletion_size: Min number of consecutive zeros in the signal for an area to be 
                           considered a gap.

        returns:
            signal_mask:         Boolean mask across signal indicating whether to include 
                                 each position in the signal.
            missing_regions:     Coordinates of signal regions with all-zero signal.
            non_missing_regions: Coordinates of signal regions with non-zero signal.
        """

        # Create array of True where signal is 0 and pad ends with False
        is_zero = np.concatenate((np.array([False]), signal == 0, np.array([False])))
        # Find coordinates of runs of zeros where signal is absent
        missing_regions = np.where(np.diff(is_zero))[0].reshape(-1, 2).astype(np.uint32)
        # Filter to find regions equal to or exceeding n base pairs
        missing_regions = missing_regions[np.where(np.diff(missing_regions).ravel() >= deletion_size)]

        # Create mask of boolean values to select signal and exclude regions of zeros
        signal_mask = np.full(signal.shape, True)

        for region in missing_regions:
            signal_mask[region[0]:region[1]] = False

        # Generate coordinates for regions with some signal
        non_missing_regions = []
        start = 0
        
        for region in missing_regions:
            if start < region[0]:
                non_missing_regions.append((start, region[0]))
            start = region[1]
        
        if start < len(signal):
            # Add end coordinates if not removed by zero gap removing filter
            non_missing_regions.append((start, len(signal)))
        
        non_missing_regions = np.array(non_missing_regions, dtype = np.uint32)

        return signal_mask, missing_regions, non_missing_regions

    def smoothSignal(self, bw_idx, chromosome, blacklist_coords = np.array([]), 
                     deletion_size = 500, kernel = np.array([]),
                     smooth_signal_file = "", missing_regions_file = "",
                     signal_stats_file = ""):
        """
        Save the original signal and smooth with a kernel.
        
        params:
            bw_idx:               ID of the sample to assess.
            chromosome:           The chromosome being analysed.
            blacklist_coords:     2D array of start and end coordinates to mask for the chromosome.
            deletion_size:        Min number of consecutive zeros in the signal for an area to be 
                                  considered a gap. 
            kernel:               Kernel for calculating convolutions across signal to create 
                                  smoothing effect. Setting kernel as None will disabled smoothing.
                                  Defaults to the recommended kernel if not set.
            smooth_signal_file:   Optional file name to save the smoothed signals.
            missing_regions_file: Optional file name to save coordinates of the missing signal / 
                                  potential deletions.
            signal_stats_file:    Set as a file to save signal statistics too. If not set, results 
                                  are returned instead.

        returns:
            signal_stats:       Dictionary of statistics for the signal. i.e. the smallest 
                                non-zero signal, number of base pairs with signal, signal sum, 
                                averages and deviations.
            smooth_signal:      Signal smoothed by the kernel.
            missing_regions:    Coordinates of regions with all-zero signal.
        """
        
        sample_names = np.array(self.getSampleNames(return_custom = True))
        sample_name = sample_names[bw_idx]
        n_samples = len(sample_names)

        # Read whole chromosome signal
        signal = self.signalReader(bw_idx = bw_idx,
                                   start_idx = 0,
                                   end_idx = -1,
                                   chromosome = chromosome,
                                   pad_end = True,
                                   verbose = self.verbose)

        if len(blacklist_coords) > 0:
            # Mask blacklisted regions
            for start, end in blacklist_coords:
                signal[start:end] = 0

        # Check whether to return results or save to file
        if smooth_signal_file or missing_regions_file or signal_stats_file:
            save_to_file = True
        else:
            save_to_file = False

        # Get signal where not zero or missing
        non_zero_mask = (signal != 0) & (~np.isnan(signal))
        # Count the number of base pairs with non-zero signal
        n_non_zero = int(np.sum(non_zero_mask))

        # Placeholder stats
        signal_stats = {"fragmentEstimate": np.nan,
                        "coverage": n_non_zero,
                        "sum": 0,
                        "mean": 0,
                        "median": 0,
                        "SD": 0,
                        "MAD": 0,
                        "meanAD": 0}

        non_zero_signal_stats = {"fragmentEstimate": np.nan,
                                 "coverage": n_non_zero,
                                 "sum": 0,
                                 "mean": np.nan,
                                 "median": np.nan,
                                 "SD": np.nan,
                                 "MAD": np.nan,
                                 "meanAD": np.nan}

        # Set placeholder array rather than large array of zeros
        smooth_signal_full = np.array([], dtype = np.float32)

        if n_non_zero == 0:
            if self.verbose > 0:
                print(f"Warning: All signal was zero for {chromosome} {sample_name}")

            if save_to_file:
                return None
            
            else:
                # No non-zero signal exists so return placeholder stats
                results = {"smooth_signal": smooth_signal_full,
                           "missing_regions": np.empty((0, 2), dtype = np.uint32),
                           "signal_stats": signal_stats,
                           "non_zero_signal_stats": non_zero_signal_stats}
                return results
            
        # Statistics for the chromosome prior to normalisation
        signal_stats = self.calculateAverages(signal = signal)
        non_zero_signal_stats = self.calculateAverages(signal = signal[non_zero_mask])

        # Calculate further stats and free memory
        non_zero_max = np.max(signal[non_zero_mask])
        non_zero_min = non_zero_signal_stats["fragmentEstimate"]

        # If some signal is positive, ensure that no negative signal exists
        mean_genome_signal = self.mean_genome_signals[bw_idx]
        contains_pos = non_zero_max > 0
        contains_neg = non_zero_min < 0
        mixed_signal = contains_neg and (mean_genome_signal > 0)
        mixed_signal = mixed_signal or (contains_pos and (mean_genome_signal < 0))

        if mixed_signal:
            print(f"Error: Cannot work with signal that is both positive and negative.\n"
                  f"{sample_name} has a minimum value of {non_zero_min} "
                  f"and a maximum value of {non_zero_max} for {chromosome}.\n" 
                  f"The mean genome signal for this sample is {mean_genome_signal}.")
            # Cannot proceed with mixed sign signal
            return None
        
        # Check if the signal is positive
        if contains_pos:
            negative = False

        # Otherwise signal must be negative
        else:
            negative = True

        if signal.dtype != "float32":
            # Convert to smaller dtype to speed up convolution
            signal = signal.astype(np.float32)

        if negative:
            # Flip the signal
            signal *= -1

        if kernel is not None:
            if len(kernel) == 0:
                # Use default kernel
                kernel = self.kernel

        if kernel is None:
            # Disable signal smoothing as user specified to not use a kernel
            apply_smoothing = False
        else:
            apply_smoothing = True
            # Centre index of kernel
            kernel_centre = int(len(kernel) / 2)

        # Find regions with consistent signal, i.e. mask consecutive runs of zeros
        signal_mask, missing_regions, non_missing_regions = self.numbaMaskMissing(signal = signal, 
                                                                                  deletion_size = deletion_size)

        # Check if any regions detected
        if non_missing_regions.shape[0] > 0:
            if (self.verbose > 0) and apply_smoothing:
                print(f"Smoothing {chromosome} signal for {sample_name} ({bw_idx + 1}/{n_samples})")

            for region_idx in range(len(non_missing_regions)):
                # Get the region coordinates, e.g. [17035200, 17035653]
                region = non_missing_regions[region_idx]
                # Set the indexes to set smoothed signal the region
                region_size = np.diff(region)[0]

                if apply_smoothing:
                    try:
                        if region_size > len(kernel):
                            # Use full kernel to smooth signal
                            signal[region[0]:region[1]] = np.convolve(signal[region[0]:region[1]],
                                                                      kernel, mode = "same")
                        else:
                            # Trim the kernel to fit the signal
                            signal[region[0]:region[1]] = np.convolve(signal[region[0]:region[1]],
                                                                      kernel[int(kernel_centre - (region_size/2)):
                                                                             int(kernel_centre + (region_size/2))],
                                                                      mode = "same")
                    except Exception as e:
                        print(f"Error smoothing {chromosome} signal:{e}")

            if apply_smoothing:
                if negative:
                    # Make smoothed signal negative again
                    signal *= -1

            if save_to_file:
                if smooth_signal_file and apply_smoothing:
                    if self.verbose > 0:
                        print(f"Saving smoothed {chromosome} signal for {sample_name} "
                              f"({bw_idx + 1}/{n_samples})") 
                    # Write smoothed signal to compressed file
                    self.saveGzip(file_name = smooth_signal_file,
                                #   array = smooth_signal_full)
                                  array = signal)

                del signal

                if missing_regions_file:
                    if self.verbose > 0:
                        print(f"Saving coordinates of missing {chromosome} signal for "
                              f"{sample_name} ({bw_idx + 1}/{n_samples})")
                    # Write missing signal regions to compressed file
                    self.saveGzip(file_name = missing_regions_file,
                                  array = missing_regions)

                del missing_regions

                if signal_stats_file:
                    if self.verbose > 0:
                        print(f"Saving {chromosome} signal statistics for {sample_name} "
                              f"({bw_idx + 1}/{n_samples})")
                    # Save signal stats to CSV
                    self.saveSignalStats(bw_idx = bw_idx,
                                         chromosome = chromosome,
                                         signal_stats = signal_stats,
                                         signal_type = "signal",
                                         file_name = signal_stats_file)
                    self.saveSignalStats(bw_idx = bw_idx,
                                         chromosome = chromosome,
                                         signal_stats = non_zero_signal_stats,
                                         signal_type = "signal_non_zero",
                                         file_name = signal_stats_file)

        else:
            del signal

            if (self.verbose > 0):
                print(f"No {chromosome} signal regions found for {sample_name} "
                      f"({bw_idx + 1}/{n_samples})")

        if not save_to_file:
            results = {#"smooth_signal": smooth_signal_full,
                       "smooth_signal": signal,
                       "missing_regions": missing_regions,
                       "signal_stats": signal_stats,
                       "non_zero_signal_stats": non_zero_signal_stats}

            return results
        
    def findCurrentFiles(self, directories, file_names_regex, single_sample = True, 
                         clear_dir_contents = False, replace_existing = False):
        r"""
        Check one or more directories to see if desired files were already created.
        
        params:
            directories:        List of directories to search,
                                e.g. ['Erythroid_Donors/Temp/Smooth_Signal/chr20', 
                                      'Erythroid_Donors/Temp/Normalised_Signals/chr20'].
            file_names_regex:   Regex to remove consistent parts of file names when extracting 
                                unique sample names,
                                e.g. "^norm-signal_*|\.npy\.gz".
                                ^ precedes the start of the file name
                                * precedes a wildcard
                                | allows the regex to match multiple cases
            single_sample:      Set as True if file names mentions only one sample, 
                                e.g. 'signals_Donor1.npy.gz'.
                                Otherwise set as False if file names mention two samples,
                                e.g. 'signals_Donor1_vs_Donor2.npy.gz'.
            clear_dir_contents: Whether to remove files in all pre-existing directories if one or 
                                more directory is missing.
            replace_existing:   Whether to overwrite previously created files.
        """

        if len(directories) == 0:
            raise ValueError("Directories cannot be empty")

        # Check to see if all directories have already been created
        directories_exist = np.array([os.path.exists(dir) for dir in directories])
        any_directories_missing = np.any(~directories_exist)

        if any_directories_missing:
            # Recreate files for all samples as any existing ones could be corrupt or from a run 
            # with different parameters
            replace_existing = True

        if not replace_existing:
            # Record pre-existing files in the output directories
            current_files = {}

        for directory, dir_idx in zip(directories, range(len(directories))):
            if directories_exist[dir_idx]:
                if clear_dir_contents and any_directories_missing:
                    # If the function is set to empty directories when only some directories exist,
                    # delete files inside the existing directory
                    for file in glob.glob(os.path.join(directory, "*")):
                        os.remove(file)

                elif not replace_existing:
                    # In the case where the directory exists but should not be cleared, record any 
                    # pre-existing files and attempt to match these with the expected files to create
                    found_files = glob.glob(os.path.join(directory, ("*")))
                    # Extract the unique sample name(s) from the files,
                    # e.g. 'Donor1' from 'signals_Donor1.npy.gz' or 
                    # 'Donor1_vs_Donor2' from 'signals_Donor1_vs_Donor2.npy.gz'
                    current_files[directory] = np.array([re.sub(file_names_regex, "", 
                                                                "_".join(file.split(os.sep)[-1].split("_")))
                                                         for file in found_files])

            else:
                # Directory not found so create empty folder to store new files
                os.makedirs(directory, exist_ok = True)

        if single_sample:
            # Set indexes for files to read signal for
            incomplete_bw_idxs = self.sample_ids
        else:
            incomplete_bw_idxs = np.arange(len(self.sample_id_pairs))

        if not replace_existing:
            # Set files to match
            if single_sample:
                target_files = self.file_sample_names
            else:
                target_files = np.array([f"{self.file_sample_names[pair[0]]}_vs_"
                                         f"{self.file_sample_names[pair[1]]}" for
                                         pair in self.sample_id_pairs])

            # Ensure names are agnostic between dashes and underscores
            target_files = [f.replace("-", "_") for f in target_files]
            # Set files to skip reading the signal for
            ignore_files = target_files.copy()
            
            for directory in current_files.keys():
                # Check to see if any of the files already exist
                check_files = [f.replace("-", "_") for f in current_files[directory]]
                ignore_files = np.intersect1d(ignore_files, check_files)

            # Set files to create as those that could not be found
            # e.g. with samples "treatment_s1.bw", "treatment_s2.bw", "control_s1.bw",
            # if files have already been created for the normalised signal as
            # "signal_treatment-s1.npy.gz", "signal_treatment-s2.npy.gz"],
            # then only a file named "signal_control-s1.npy.gz" needs to be created
            incomplete_bw_idxs = np.setdiff1d(incomplete_bw_idxs,
                                              np.isin(target_files, ignore_files).nonzero()[0])

        return incomplete_bw_idxs

    def interleaveChroms(self, process_chrom_bws):
        """
        Create a list of (chromosome, sample_id) pairs to optimise process running order. 
        Chromosomes are interleaved by size, so that the remaining largest is followed by 
        the remaining smallest, which is then followed by the next remaining largest etc.
        e.g. [("chr1", 1), ("chrY", 1"), ("chr1", 2), ("chrY", 2), ("chr2", 1), ("chr21", 1), 
              ("chr2", 2), ("chr21", 2)].

        params:
            process_chrom_bws: List of pyBigWig bigWigs to process.

        returns:
            interleaved_chroms: List of (chromosome, sample_id) pairs in order to process.
        """

        # Get chromosomes from bigWigs to process and sort by size
        sorted_chroms = sorted(process_chrom_bws,
                               key = lambda x: self.chrom_sizes.get(x[0], float("inf")))

        # Interleave largest and smallest chromosomes to reduce memory usage
        interleaved_chroms = []
        left = 0
        right = len(sorted_chroms) - 1
        use_largest = True

        while left <= right:
            if use_largest:
                # Add the largest remaining chromosome-sample tuple
                interleaved_chroms.append(sorted_chroms[right])
                right -= 1
            else:
                # Add the smallest remaining chromosome-sample tuple
                interleaved_chroms.append(sorted_chroms[left])
                left += 1

            # Flip order
            use_largest = not use_largest

        return interleaved_chroms

    def plotKernel(self, kernel = [], title = ""):
        """
        Create a plot of a smoothing kernel.

        params:
            kernel: List or array containing kernel to plot. By default, the set kernel is shown.
            title:  Set a custom title for the plot.
        """

        if len(kernel) == 0:
            kernel = self.kernel
        else:
            kernel = np.array(kernel)

        if len(title) == 0:
            title = f"{self.analysis_name.replace('_', ' ')} Smoothing Kernel"

        plt.figure(figsize = (4, 3))
        plt.plot(kernel)
        plt.title(title, fontweight = "bold")
        plt.xlabel("Position (Base Pairs)")
        plt.ylabel("Weight")
        plt.show()

    def convolveSignals(self, chromosomes = [], replace_existing = False):
        """
        Perform initial steps prior to signal zone prediction. These include reading signal 
        from each bigWig, applying convolution to create a smoothed version of this signal,
        calculating summary statistics about the signal, identifying potential deletions and 
        removal of blacklisted regions.

        params:
            chromosomes:      List of chromosomes to read signal from per bigWig.
            replace_existing: Whether to overwrite previously created files.
        """

        if len(chromosomes) == 0:
            chromosomes = self.chromosomes
        else:
            missing_chroms = [chrom for chrom in chromosomes if chrom not in self.chromosomes]

            if missing_chroms:
                raise ValueError(f"{len(missing_chroms)} chromosomes were not found: "
                                 f"{', '.join(missing_chroms)}")

        # Create tuples of chromosomes and samples to process
        process_chrom_bws = []

        for chrom in chromosomes:
            chrom_smooth_signal_dir = os.path.join(self.output_directories["smooth_signal"], chrom)
            chrom_sample_missing_dir = os.path.join(self.output_directories["missing_signal"], chrom)
            chrom_signal_stats_dir = os.path.join(self.output_directories["signal_stats"], chrom)

            # Find indexes of samples where signals, smoothed signal and deletions have not yet 
            # been saved to file
            incomplete_bw_idxs = self.findCurrentFiles(directories = [chrom_smooth_signal_dir,
                                                                      chrom_sample_missing_dir,
                                                                      chrom_signal_stats_dir],
                                                       file_names_regex = "^smooth-signal_*|" +
                                                                          "^missing-regions_*|" + 
                                                                          r"\.npy\.gz|" +
                                                                          "^signal-stats_*|" + 
                                                                          r"\.csv",
                                                       replace_existing = replace_existing)
            
            for bw_idx in incomplete_bw_idxs:
                process_chrom_bws.append((chrom, bw_idx))

        if (self.n_cores > 1) and self.interleave_sizes:
            process_chrom_bws = self.interleaveChroms(process_chrom_bws)

        if len(process_chrom_bws) > 0:
            # Replace any underscores or dots as these are excluded from file names
            file_sample_names = self.file_sample_names

            # Blacklist coordinate arrays per chromosome
            blacklist_coords = {}

            for chrom in chromosomes:
                if self.blacklist is not None:
                    blacklist_coords[chrom] = np.array(self.getBlacklist(chromosome = chrom)[["start", "end"]])
                else:
                    blacklist_coords[chrom] = np.array([])

            if self.n_cores > 1:
                with ProcessPoolExecutor(self.n_cores) as executor:
                    smooth_processes = []
                    for chrom, bw_idx in process_chrom_bws:
                        # Smooth signal and fit distributions to identify regions of the chromosome with signal
                        smooth_processes.append(executor.submit(self.smoothSignal,
                                                                bw_idx = bw_idx,
                                                                chromosome = chrom,
                                                                blacklist_coords = blacklist_coords[chrom],
                                                                deletion_size = self.deletion_size,
                                                                smooth_signal_file = os.path.join(self.output_directories["smooth_signal"], chrom,
                                                                                                  f"smooth-signal_{file_sample_names[bw_idx]}"),
                                                                missing_regions_file = os.path.join(self.output_directories["missing_signal"], chrom,
                                                                                                    f"missing-regions_{file_sample_names[bw_idx]}"),
                                                                signal_stats_file = os.path.join(self.output_directories["signal_stats"], chrom,
                                                                                                 f"signal-stats_{file_sample_names[bw_idx]}.csv")))
                    
                    if self.checkParallelErrors(smooth_processes):
                        return None

            else:
                for chrom, bw_idx in process_chrom_bws:
                    self.smoothSignal(bw_idx = bw_idx,
                                      chromosome = chrom,
                                      blacklist_coords = blacklist_coords[chrom],
                                      deletion_size = self.deletion_size,
                                      smooth_signal_file = os.path.join(self.output_directories["smooth_signal"], chrom,
                                                                        f"smooth-signal_{file_sample_names[bw_idx]}"),
                                      missing_regions_file = os.path.join(self.output_directories["missing_signal"], chrom,
                                                                          f"missing-regions_{file_sample_names[bw_idx]}"),
                                      signal_stats_file = os.path.join(self.output_directories["signal_stats"], chrom,
                                                                       f"signal-stats_{file_sample_names[bw_idx]}.csv"))

        elif self.verbose > 0:
            print(f"Signal, smoothed signal and missing signal already created for chromosomes")

    def readArrayFile(self, bw_idx, directory, file_prefix):
        """
        Read and return a numpy array from a file containing a sample name.
        
        params:
            bw_idx:      Index of the sample the file is linked to.
            directory:   Directory to read file from.
            file_prefix: Shared beginning of file name, e.g. "signal" in "signal_1.npy.gz".
        
        returns:
            bw_idx: Index of the sample the file is linked to (used as an identifier 
                    during parallelisation).
            data:   The numpy array read from file.
        """

        if file_prefix:
            file_prefix += "_"

        sample_file = os.path.join(directory, (file_prefix + self.file_sample_names[bw_idx] + ".npy.gz"))
        if os.path.isfile(sample_file):
            if os.path.getsize(sample_file) > 0:
                data = self.loadGzip(file_name = sample_file)
            else:
                raise ValueError(f'File "{sample_file}" is empty')
        else:
            raise FileNotFoundError(f"Cannot find file {sample_file}")

        return bw_idx, data

    def readSignalFile(self, bw_idx, chromosome, signal_zones, read_normalised = True,
                       zone_type = None, verbose = 0, dtype = object):
        """
        Read regions of a chromosome signal from file.
        
        params:
            bw_idx:          Index of the sample the file is linked to.
            chromosome:      Chromosome to read signal from, e.g. "chr1".
            signal_zones:    List of coordinate pairs for regions to keep signal for.
            read_normalised: Whether to read the original signal for the sample (set as False) or 
                             the normalised signal (set as True).
            zone_type:       Set as either 'padded' or 'unpadded' to match the type of signal_zones. 
                             Not required, but can be used to print progress.
            verbose:         Set as a number > 0 to print progress.
            dtype:           Datatype to return results as.
        """

        # Get the name of the sample according to its ID
        custom_sample_name = self.getSampleNames(return_custom = True)[bw_idx]
        file_sample_name = self.file_sample_names[bw_idx]

        if verbose > 0:
            print(f"Extracting {self.norm_method + ' normalised' if read_normalised else 'original'} "
                  f'{zone_type + " " if zone_type is not None else ""}zones of signal ' 
                  f"for {chromosome} {custom_sample_name} ({bw_idx + 1}/{len(self.sample_names)})")

        if read_normalised:
            try:
                # Open normalised signal from file
                signal = self.loadGzip(file_name = os.path.join(self.output_directories["norm_signal"], 
                                                                self.norm_method, chromosome,
                                                                f"{self.norm_method}-norm-signal_{file_sample_name}.npy.gz"))
            except Exception as e:
                print(f"Could not open {self.norm_method} normalised signal for {custom_sample_name} "
                      f"due to exception:\n{e}")
                return None, None
        else:
            try:
                # Read whole chromosome signal
                signal = self.signalReader(bw_idx = bw_idx,
                                           start_idx = 0,
                                           end_idx = -1,
                                           pad_end = True,
                                           chromosome = chromosome)
            except Exception as e:
                print(f"Could not open original signal for {custom_sample_name} due to exception:\n{e}")
                return None, None
                
        # Filter to get signal within ranges
        signal = np.array([signal[coords[0]:coords[1]] for coords in signal_zones], dtype = dtype)
        
        return bw_idx, signal

    def saveGzip(self, file_name, array):
        """ 
        Write numpy array to gzipped file.
        
        params:
            file_name: Name of file to save results to.
            array:     Numpy array of results to save.
        """

        if not file_name.endswith(".npy.gz"):
            # Add numpy and gzip file extension
            file_name += ".npy.gz"

        # Save array to a compressed file
        with gzip.GzipFile(file_name, "wb") as gzip_file:
            np.save(file = gzip_file, arr = array)

    def loadGzip(self, file_name):
        """
        Open numpy array from gzipped file.

        params:
            file_name: Name of file to open results from.

        returns:
            data: Contents of the file.
        """

        if not file_name.endswith(".npy.gz"):
            # Add numpy and gzip file extension
            file_name += ".npy.gz"

        # Open array in a compressed file
        with gzip.GzipFile(file_name, "rb") as gzip_file:
            data = np.load(file = gzip_file, allow_pickle = True)

        return data

    def testDistributions(self, chromosomes = [], log_transform = None, downsample_size = None, 
                          replace_existing = False):
        """
        Test which distributions fit the data best to use for signal zone prediction.

        params:
            chromosomes:      List of chromosomes.
            log_transform:    Whether to apply log transformation to signal before fitting 
                              distribution(s).
            downsample_size:  Number of positions of signal to select.
            replace_existing: Whether to overwrite previously created files.
        """

        if len(chromosomes) == 0:
            chromosomes = self.chromosomes

        # Create tuples of chromosomes and samples to process
        process_chrom_bws = []

        for chrom in chromosomes:
            signal_stats_folder = os.path.join(self.output_directories["signal_stats"], chrom)

            if not os.path.exists(signal_stats_folder):
                # Cannot test distributions if signal statistics not yet created
                raise FileNotFoundError(f"No signal statistics yet calculated for {chrom}.\n"
                                        f"First run convolveSignals.")
            
            chrom_dist_stats_dir = os.path.join(self.output_directories["dist_stats"], chrom)

            # Find indexes of samples where distribution statistics are missing
            incomplete_bw_idxs = self.findCurrentFiles(directories = [chrom_dist_stats_dir],
                                                       file_names_regex = "^distribution-stats_*|" + 
                                                                          r"\.csv",
                                                       replace_existing = replace_existing)
            
            for bw_idx in incomplete_bw_idxs:
                process_chrom_bws.append((chrom, bw_idx))

        if (self.n_cores > 1) & self.interleave_sizes:
            process_chrom_bws = self.interleaveChroms(process_chrom_bws)

        if len(process_chrom_bws) > 0:
            # Update attributes if different parameters given
            if downsample_size is None:
                downsample_size = self.downsample_size
            # No downsampling if set as zero
            elif downsample_size != 0:
                self.setDownsampleSize(downsample_size, verbose = self.verbose)

            if log_transform is None:
                log_transform = self.log_transform
            else:
                self.setLogTransform(log_transform, verbose = self.verbose)

            if (self.n_cores > 1) and (len(self.sample_ids) > 1):
                # Predict regions of signal per sample from the smoothed signal
                with ProcessPoolExecutor(self.n_cores) as executor:
                    test_processes = []
                    for chrom, bw_idx in process_chrom_bws:
                        # Set files
                        file_sample_name = self.file_sample_names[bw_idx]
                        signal_stats_file = os.path.join(signal_stats_folder,
                                                         f"signal-stats_{file_sample_name}.csv")
                        dist_stats_file = os.path.join(self.output_directories["dist_stats"], chrom,
                                                       f"distribution-stats_{file_sample_name}.csv")
                        
                        # Run process
                        test_processes.append(executor.submit(self.testZoneCallDistributions,
                                                              bw_idx = bw_idx,
                                                              chromosome = chrom,
                                                              log_transform = log_transform,
                                                              downsample_size = downsample_size,
                                                              signal_stats_file = signal_stats_file,
                                                              dist_stats_file = dist_stats_file))

                    # Check for errors during multiprocessing
                    self.checkParallelErrors(test_processes)

            else:
                # Run sequentially rather than using parallelisation
                for chrom, bw_idx in process_chrom_bws:
                    file_sample_name = self.file_sample_names[bw_idx]
                    signal_stats_file = os.path.join(signal_stats_folder,
                                                     f"signal-stats_{file_sample_name}.csv")
                    dist_stats_file = os.path.join(self.output_directories["dist_stats"], chrom,
                                                   f"distribution-stats_{file_sample_name}.csv")
                    self.testZoneCallDistributions(bw_idx = bw_idx,
                                                   chromosome = chrom,
                                                   log_transform = log_transform,
                                                   downsample_size = downsample_size,
                                                   signal_stats_file = signal_stats_file,
                                                   dist_stats_file = dist_stats_file)

        elif self.verbose > 0:
            print(f"Distribution statistics already calculated for chromosomes")

    def calculateAverages(self, signal):
        """
        Create dictionary of average and deviation statistics for a signal.

        params:
            signal: Array of signal to calculate statistics for.
        """

        stats = {"fragmentEstimate": np.nan,
                 "coverage": np.nan}

        if len(signal) == 0:
            # If signal is missing, return placeholder stats
            stats = {"fragmentEstimate": np.nan,
                     "coverage": np.nan,
                     "sum": 0,
                     "mean": np.nan,
                     "median": np.nan,
                     "MAD": np.nan,
                     "meanAD": np.nan,
                     "SD": np.nan}

            return stats

        non_zero_mask = (signal != 0) & (~np.isnan(signal))

        if np.sum(non_zero_mask) > 0:
            min_non_zero = min(signal[non_zero_mask])
            max_non_zero = max(signal[non_zero_mask])

            signal_sign = "mixed"

            if min_non_zero > 0:
                signal_sign = "positive"
            elif max_non_zero < 0:
                signal_sign = "negative"

            # These statstics can only be calculated when the signal has the same sign
            if signal_sign != "mixed":
                if signal_sign == "negative":
                    stats["fragmentEstimate"] = float(max(signal[non_zero_mask]))
                else:
                    stats["fragmentEstimate"] = float(min(signal[non_zero_mask]))

                # Number of base pairs of non-zero signal
                stats["coverage"] = int(np.sum(non_zero_mask))

        # Absolute sum
        stats["sum"] = float(np.nansum(np.abs(signal)))
        # Averages and deviations
        stats["mean"] = float(np.nanmean(signal))
        stats["median"] = float(np.nanmedian(signal))
        stats["MAD"] = float(np.nanmedian(np.abs(signal - stats["median"])))
        stats["meanAD"] = float(np.nanmean(np.abs(signal - stats["mean"])))
        stats["SD"] = float(np.nanstd(signal))

        return stats

    def transformSignal(self, signal, log_transform = True, downsample_size = 0, 
                        calculate_stats = False):
        """
        Transform the smoothed / original signal before fitting distibutions for zone prediction.

        params:
            signal:          Array of signal to transform.
            log_transform:   Whether to apply log transformation to signal.
            downsample_size: Set as an integer > 0 to downsample signal by this many base pairs.
            calculate_stats: Whether to calculate statistics from the transformed signal.

        returns:
            signal_to_fit:     The transformed signal.
            transformed_stats: Dictionary of statistics from the transformed signal. 
                               Only returned if calculate_stats is True.
        """

        # Mask zeros and log transform to reduce right-skew and zero-inflation
        non_zero_mask = np.where((signal != 0) & (~np.isnan(signal)))
        # Get absolute value to prevent signal being negative
        non_zero_signal = np.abs(signal[non_zero_mask])
        del signal

        if log_transform:
            transformed_signal = np.log(non_zero_signal)
        else:
            transformed_signal = non_zero_signal

        if (downsample_size > 0) and (downsample_size < len(transformed_signal)):
            downsample = True
        else:
            downsample = False

        # Calculate statistics for the transformed signal
        if calculate_stats:
            # Create dictionary of statistics
            transformed_stats = self.calculateAverages(transformed_signal)

            # Set statistics that could not be set due to the log transformation
            if log_transform:
                # Calculate estimated value for a single read
                fragment_estimate = min(non_zero_signal)
                fragment_estimate = np.log(fragment_estimate)
                transformed_stats["fragmentEstimate"] = fragment_estimate
                # Count number of non-zero base pairs
                transformed_stats["coverage"] = len(non_zero_signal)

            if downsample:
                # Instead use number of downsampled base pairs
                transformed_stats["coverage"] = downsample_size

        del non_zero_signal

        if downsample:
            # Seed for reproducibility
            np.random.seed(self.downsample_seed)
            # Subsample signal for speed up in distribution fitting
            signal_to_fit = np.random.choice(transformed_signal, size = downsample_size, 
                                             replace = False)
        else:
            signal_to_fit = transformed_signal

        if calculate_stats:
            return signal_to_fit, transformed_stats
        else:
            return signal_to_fit

    def calculateMeanScale(self, dist_name, sd):
        """
        Adjust standard deviation to scale that approximates SD for normal distribution.

        params:
            dist_name: Name of distribution to fit.
            sd:        Value of one standard deviation.
        """

        # Set the scale factor (default no adjustment)
        k = 1

        if dist_name == "laplace":
            # Laplace scale is b
            # sd = b * sqrt(2) -> k = 1 / sqrt(2)
            k = 1 / np.sqrt(2)

        elif dist_name == "logistic":
            # Logistic scale is s
            # sd = (s * pi) / sqrt(3) -> k = sqrt(3) / pi
            k = np.sqrt(3) / np.pi

        elif dist_name == "gumbel_l" or dist_name == "gumbel_r":
            # Gumbel scale is beta
            # sd = (pi / sqrt(6)) * beta -> k = sqrt(6) / pi
            k = np.sqrt(6) / np.pi

        # Set scale as standard deviation adjusted by k
        scale = k * sd

        return scale

    def calculateMedianScale(self, dist_name, mad, meanAD):
        """
        Translate MAD or mean average deviation to scale that approximates SD for normal distribution.

        params:
            dist_name: Name of distribution to fit.
            mad:       Value of one median absolute deviation.
            meanAD:    Value of one mean average deviation.
        """
        
        # Set the scale factor (default no adjustment)
        k = 1

        if mad == 0:
            if dist_name == "norm":
                # Normal scale is sigma, meanAD derived from the mean of the half-normal distribution
                # MeanAD = sigma * sqrt(2/pi) -> k = 1 / sqrt(2 / pi)
                k = 1 / np.sqrt(2 / np.pi)

            elif dist_name == "logistic":
                # Logistic scale is s
                # MeanAD = 2s * ln(2) -> k = 1 / 2ln(2)
                k = 1 / (2 * np.log(2))

            elif dist_name == "gumbel_l" or dist_name == "gumbel_r":
                # Gumbel scale is beta
                # MeanAD = G - Ei -> k ≈ 1.6654

                # Get Euler's number and Euler-Mascheron constant
                euler_e = sp.E
                euler_gamma = np.euler_gamma
                # Meijer G-function for G^{2,0}_{2,1}(e, gamma^-1 | (1, 1), (0))
                z = euler_e / euler_gamma
                meijer_g = sp.functions.special.hyper.meijerg([1,1], [], [], [0], z).evalf()
                # Exponential integral for Ei(-e^(-gamma))
                ei_value = expi(-np.exp(-euler_gamma))
                # Compute the constant from the difference
                k = float(meijer_g - ei_value)

            # Set scale as mean absolute deviation adjusted by k
            scale = k * meanAD

        else:
            # Use MAD as approximation
            if dist_name == "norm":
                # Normal scale is sigma, phi^-1 is the inverse of the CDF of the normal distribution
                # MAD = sigma * phi^-1(3/4) -> k ≈ 1.4826
                k = 1.4826

            elif dist_name == "laplace":
                # Laplace scale is b
                # MAD = b * ln(2) -> k = 1 / ln(2)
                k = 1 / np.log(2)

            elif dist_name == "gumbel_l" or dist_name == "gumbel_r":
                # Gumbel scale is beta and 0.7670 dervied from numerical approximation
                # MAD ≈ beta * 0.7670 -> k ≈ 1.3036
                k = 1 / 0.767049251325708

            # Set scale as median absolute deviation adjusted by k
            scale = k * mad

        return scale

    def testZoneCallDistributions(self, bw_idx, chromosome, log_transform = True, downsample_size = 300,
                                  signal_stats_file = "", dist_stats_file = ""):
        """
        Tests fitted distributions on smoothed signal to later determine which to use to predict 
        signal zones.
        
        params:
            bw_idx:               Sample ID to test distributions for.
            chromosome:           Chromosome to test.
            log_transform:        Whether to log transform the signal before distribution fitting.
            downsample_size:      Number of positions of signal to select.
            signal_stats_file:    Set as a file path to calculate signal statistics and save these 
                                  to file.
            dist_stats_file:      Set as a file path to calculate distribution statistics and save 
                                  these to file.
        """
        
        sample_name = np.array(self.getSampleNames(return_custom = True))[bw_idx]
        n_samples = len(self.sample_ids)

        # Check whether to return results or save to file
        if signal_stats_file or dist_stats_file:
            save_to_file = True
        else:
            save_to_file = False

        # Test if smoothing was applied
        if self.kernel is not None:
            smoothed = True
            try:
                # Open smoothed signal for the sample's chromosome
                signal = self.getSmoothedSignal(sample = bw_idx, chromosome = chromosome)

            except Exception as e:
                print(f"Warning: Could not read smoothed signal for {chromosome} {sample_name} "
                      f"due to exception\n")
                raise e

        else:
            smoothed = False
            try:
                # Open full signal prior to normalisation
                signal = self.signalReader(bw_idx = bw_idx,
                                           start_idx = 0,
                                           end_idx = -1,
                                           chromosome = chromosome,
                                           pad_end = True,
                                           verbose = self.verbose)
            except:
                # If signal is missing, analysis may not have been run on the sample's chromosome
                raise ValueError(f"Could not read signal for {chromosome} {sample_name} due to "
                                 f"exception\n{e}")

        if (self.verbose > 0):
            print(f'Transforming {"smoothed" if smoothed else ""} signal',
                  f"for {chromosome} {sample_name} ({bw_idx + 1}/{n_samples})")

        # Transform the signal before distribution fitting and calculate stats before downsampling
        signal_to_fit, transformed_stats = self.transformSignal(signal = signal,
                                                                log_transform = log_transform,
                                                                downsample_size = downsample_size,
                                                                calculate_stats = True)

        del signal
        gc.collect()

        if signal_stats_file:
            if self.verbose > 0:
                print(f'Saving transformed {"smoothed" if smoothed else ""} signal',
                      f"statistics for {chromosome} {sample_name} ({bw_idx + 1}/{n_samples})")
            # Save transformed signal stats to CSV
            self.saveSignalStats(bw_idx = bw_idx,
                                 chromosome = chromosome,
                                 signal_stats = transformed_stats,
                                 signal_type = "signal_transformed",
                                 file_name = signal_stats_file)

        if (self.verbose > 0):
            print(f"Fitting distributions to {chromosome} {sample_name} ({bw_idx + 1}/{n_samples})")

        if downsample_size > 0:
            # Update to get statistics for the downsampled statistics
            transformed_stats = self.calculateAverages(signal_to_fit)

        # Statistics for all fitted distributions
        distribution_stats = {}

        for dist_name in self.test_distributions:
            # Get the distribution class
            dist = self.scipy_distributions[dist_name]

            # Stats for the individual distribution
            dist_stats = {}

            # Test different combinations of parameters
            for param_type in self.param_types:
                try:
                    if param_type == "scipy_fit":
                        # Fit distribution to data to estimate parameters
                        params = dist.fit(signal_to_fit)

                    else:
                        if param_type == "mean_fit":
                            # Set estimated location
                            location = transformed_stats["mean"]
                            # Set estimated scale
                            scale = self.calculateMeanScale(dist_name = dist_name,
                                                            sd = transformed_stats["SD"])
                        else:
                            location = transformed_stats["median"]
                            scale = self.calculateMedianScale(dist_name = dist_name,
                                                              mad = transformed_stats["MAD"],
                                                              meanAD = transformed_stats["meanAD"])
                                    
                        # Set manual parameter estimates
                        params = (location, scale)

                    # Apply Kolmogorov-Smirnov goodness of fit test
                    ks_stat, ks_p_value = stats.kstest(signal_to_fit, dist.cdf, args = params)
                        
                    # Compute log-likelihood, AIC and BIC
                    signal_prob = dist.pdf(signal_to_fit, *params)

                    if np.any(np.isin(signal_prob, 0)):
                        log_likelihood = - np.inf
                        aic = np.inf
                        bic = np.inf
                    else:
                        log_likelihood = np.sum(np.log(signal_prob))
                        n_params = len(params)
                        aic = 2 * n_params - 2 * log_likelihood
                        bic = n_params * np.log(len(signal_to_fit)) - 2 * log_likelihood
                    
                    dist_stats[param_type] = {"location": params[0],
                                              "scale": params[1],
                                              "KS_stat": ks_stat,
                                              "KS_p_value": ks_p_value,
                                              "log_likelihood": log_likelihood,
                                              "AIC": aic,
                                              "BIC": bic}

                except Exception as e:
                    if self.verbose > 0:
                        print(f"Error fitting {dist_name} {param_type} parameters for {chromosome} "
                              f"{sample_name}")
                        print(e)
                
            if len(dist_stats) > 0:
                # Record only if the distribution was fitted
                distribution_stats[dist_name] = dist_stats

        if save_to_file:
            if dist_stats_file:
                if self.verbose > 0:
                    print(f"Saving distribution fitting statistics for {chromosome} {sample_name} "
                          f"({bw_idx + 1}/{n_samples})")
                # Save distribution stats to CSV
                self.saveDistributionStats(distribution_stats = distribution_stats,
                                           file_name = dist_stats_file)
        else:
            # Return dictionary of statistics
            return transformed_stats, distribution_stats
        
    def testBestByStat(self, stats_df, stat_column, name_column):
        """
        Find the top one or two performing entrys in a DataFrame according to a specified statistic.

        params:
            stats_df:    DataFrame with statistics and labels.
            stat_column: Column containing statistic to sort by, e.g. "KS_stat".
            name_column: Column containing label of entry for statistic, e.g. "distribution" column.
        """

        if stat_column in ["KS_p_value"]:
            # Sort highest to lowest
            ascending = False
        else:
            # Sort lowest to highest
            ascending = True

        # Find the combinations that scored best for the statistic per chromosome
        sorted_stats_df = stats_df.sort_values(by = stat_column, ascending = ascending)
        sorted_stats_df = sorted_stats_df.drop_duplicates(subset = "chrom", keep = "first")

        if sorted_stats_df[name_column].nunique() > 1:
            # Identify the top two types according to the statistic
            column_count = sorted_stats_df[name_column].value_counts()
            count_1 = column_count.iloc[0]
            count_2 = column_count.iloc[1]
            name_1 = column_count.index[0]
            name_2 = column_count.index[1]
        else:
            # Only one type found
            name_1 = np.unique(sorted_stats_df[name_column])[0]
            name_2 = ""
            count_1 = len(sorted_stats_df)
            count_2 = 0

        return name_1, name_2, count_1, count_2

    def findBestColumnValue(self, stats_df, name_column, primary_stat, secondary_stat = "", 
                            tertiary_stat = ""):
        """
        Find the best result within a DataFrame first by a primary statistic, or if tied by a 
        secondary or tertiary statistic.

        params:
            stats_df:       DataFrame with one to three statistics and labels.
            name_column:    Column containing label of entry for statistics, e.g. "distribution" 
                            column.
            primary_stat:   Column containing first statistic to sort by.
            secondary_stat: Column containing optional second statistic to sort by.
            tertiary_stat:  Column containing optional third statistic to sort by.
        """

        # Find the combinations that scored best for the primary statistic per chromosome
        name_1, name_2, count_1, count_2 = self.testBestByStat(stats_df,
                                                               primary_stat,
                                                               name_column)

        # If tied, repeat for the secondary and tertiary statistics
        if secondary_stat != "" and count_1 == count_2:
            name_1, name_2, count_1, count_2 = self.testBestByStat(stats_df,
                                                                   secondary_stat,
                                                                   name_column)
            if tertiary_stat != "" and count_1 == count_2:
                name_1, name_2, count_1, count_2 = self.testBestByStat(stats_df,
                                                                       tertiary_stat,
                                                                       name_column)

        if count_1 == count_2:
            # Record both as best due to a repeated tie
            best_name = [name_1, name_2]
        else:
            # Counts were ordered, so select the name which had the most
            best_name = [name_1]

        return best_name

    def inferBestDistribution(self, dist_stats_df = None, zone_distribution = None, 
                              primary_stat = "KS_stat", secondary_stat = "AIC", 
                              tertiary_stat = "BIC", use_rank = False):
        """
        Estimate which distribution and parameter type best fit across all the samples.

        params:
            dist_stats_df:  DataFrame of combine distribution statistics across all samples.
            primary_stat:   Column containing first statistic to sort by.
            secondary_stat: Column containing optional second statistic to sort by.
            tertiary_stat:  Column containing optional third statistic to sort by.
            use_rank:       If True ranking is used to infer the best. Otherwise, the average of 
                            the primary statistic is used.
        """

        custom_sample_names = np.array(self.getSampleNames(return_custom = True))

        if dist_stats_df is None:
            dist_stats_file = os.path.join(self.output_directories["output_stats"], 
                                           "combined-distribution-stats.csv")

            if not os.path.isfile(dist_stats_file):
                # Create combined distribution statistics if does not yet exist
                dist_stats_df = self.combineStats(stats_type = "distribution",
                                                  file_name = dist_stats_file,
                                                  return_df = True)
            else:
                # Read from file
                dist_stats_df = pd.read_csv(dist_stats_file)

        elif not isinstance(dist_stats_df, pd.DataFrame):
            raise ValueError("dist_stats_df must be given as a DataFrame or set to None to"
                            "generate combined distribution statistics")

        if zone_distribution is not None:
            # Check that distribution exists
            if not ((dist_stats_df["distribution"] == zone_distribution).any()):
                raise ValueError(f"Could not find distribution {zone_distribution} in"
                                 f"distribution statistics")

            # Disable writing best parameters as distribution may not be best overall fit
            store_best_attribute = False
            # Filter for rows with the distribution
            dist_stats_df = dist_stats_df[dist_stats_df["distribution"] == zone_distribution]
        else:
            store_best_attribute = True
            # Filter for rows with set distributions to test
            dist_stats_df = dist_stats_df[np.isin(dist_stats_df["distribution"], self.test_distributions)]

        # Count how many distributions were fitted and the types of parameter combinations
        n_distributions = dist_stats_df["distribution"].nunique()
        n_param_types = dist_stats_df["param_type"].nunique()

        # Check if only one distribution was fitted with one parameter type
        if (n_distributions == 1) and (n_param_types == 1):
            best_distribution = dist_stats_df["distribution"].unique()[0]
            best_param_type = dist_stats_df["param_type"].unique()[0]

            if self.verbose > 0:
                print(f"Only the {best_distribution} distribution was fitted with {best_param_type} "
                      f"parameters")

            # Convert to dictionary
            best_params = {best_distribution: [best_param_type]}

            if store_best_attribute:
                # Set attribute to the distribution and parameter type
                self.best_zone_distribution = best_params

            return best_params

        # Ensure the specified statistics exist in the columns
        for stat in [primary_stat, secondary_stat, tertiary_stat]:
            if stat not in dist_stats_df.columns:
                raise ValueError(f'Error: "{stat}" not found in DataFrame columns: '
                                 f"{list(dist_stats_df.columns)}")

        # Check if any samples are not in the distribution statistics
        missing_samples = np.setdiff1d(custom_sample_names, dist_stats_df["sample"].values)

        if len(missing_samples) > 0:
            raise ValueError(f"The following {len(missing_samples)} samples are missing in "
                             f"the distribution statistics: {','.join(missing_samples)}")

        # Columns in distribution statistics to calculate averages for
        results_columns = list(set(["AIC", "BIC", "KS_p_value", "KS_stat", "log_likelihood"]) & 
                               set(dist_stats_df.columns))

        if primary_stat in ["KS_p_value"]:
            # Sort highest to lowest
            ascending = False
        else:
            # Sort lowest to highest
            ascending = True

        # Calculate the mean statistics per sample for each distribution and parameter type
        sample_averages = dist_stats_df.groupby(["sample",
                                                 "distribution",
                                                 "param_type"])[results_columns].mean().reset_index()
        # Average across samples
        dist_averages = sample_averages.groupby(["distribution",
                                                 "param_type"])[results_columns].mean().reset_index()
        dist_averages = dist_averages.sort_values(by = primary_stat, ascending = ascending)

        if store_best_attribute:
            # Save averaged distribution statistics to file
            mean_dist_stats_file = os.path.join(self.output_directories["output_stats"], 
                                                "mean-distribution-stats.csv")
            dist_averages.to_csv(mean_dist_stats_file, header = True, index = False)

        if use_rank:
            if n_distributions > 1:
                # Find the best distribution per sample
                best_sample_dists = {}
                all_best_dist = []

                for bw_idx in self.sample_ids:
                    # Get rows for the sample
                    custom_sample_name = custom_sample_names[bw_idx]
                    sample_dist_df = dist_stats_df.loc[dist_stats_df["sample"] == custom_sample_name]
                    # Find the name of the distribution that performs best the most times
                    dist_name = self.findBestColumnValue(stats_df = sample_dist_df,
                                                         name_column = "distribution",
                                                         primary_stat = primary_stat,
                                                         secondary_stat = secondary_stat,
                                                         tertiary_stat = tertiary_stat)
                    all_best_dist.append(dist_name)
                    best_sample_dists[custom_sample_name] = dist_name

                if len(np.unique(all_best_dist)) == 1:
                    # One distribution consistently performed best
                    best_distribution = all_best_dist[0]
                else:
                    # Otherwise find the one that was best across samples
                    column_count = pd.Series(all_best_dist).value_counts()
                    count_1 = column_count.iloc[0]
                    count_2 = column_count.iloc[1]
                    name_1 = column_count.index[0]

                    if count_1 == count_2:
                        # Tied so record both
                        name_2 = column_count.index[1]
                        best_distribution = [name_1, name_2]
                    else:
                        best_distribution = name_1

            else:
                # Only one distribution fitted so set it as best for all samples and best overall
                best_distribution = list(np.unique(dist_stats_df["distribution"]))

            if n_param_types > 1:
                best_sample_params = {}
                all_best_params = dict(zip(best_distribution, [[]] * len(best_distribution)))
                best_params = {}

                # Find the best parameter type per sample
                for bw_idx in self.sample_ids:
                    custom_sample_name = custom_sample_names[bw_idx]
                    best_sample_params[custom_sample_name] = {}
                    # Get rows for the sample
                    sample_dist_df = dist_stats_df.loc[dist_stats_df["sample"] == custom_sample_name]

                    for dist_name in best_distribution:
                        # Find the parameter type that performs best the most times
                        param_type = self.findBestColumnValue(stats_df = sample_dist_df,
                                                              name_column = "param_type",
                                                              primary_stat = primary_stat,
                                                              secondary_stat = secondary_stat,
                                                              tertiary_stat = tertiary_stat)
                        all_best_params[dist_name].append(param_type)
                        best_sample_params[custom_sample_name][dist_name] = param_type

                for dist_name in best_distribution:
                    if len(np.unique(all_best_params[dist_name])) == 1:
                        # One parameter type consistently performed best
                        best_params[dist_name] = all_best_params[dist_name][0]
                    else:
                        # Otherwise find the one that was best across samples
                        column_count = pd.Series(all_best_params[dist_name]).value_counts()
                        count_1 = column_count.iloc[0]
                        count_2 = column_count.iloc[1]
                        name_1 = column_count.index[0]
                        if count_1 == count_2:
                            # Tied so record both
                            name_2 = column_count.index[1]
                            best_params[dist_name] = [name_1[0], name_2[0]]
                        else:
                            best_params[dist_name] = name_1

            else:
                # Only one parameter type tested
                best_param_type = [dist_stats_df["param_type"].loc[0]]
                best_params = dict(zip(best_distribution, [best_param_type] * len(best_distribution)))

        else:
            # Get row with the best distribution and parameter type
            top_row = dist_averages.iloc[0]
            best_params = {top_row["distribution"]: [top_row["param_type"]]}

        if store_best_attribute:
            # Set inferred attribute as best distribution(s) and parameter types
            self.best_zone_distribution = best_params

        return best_params

    def calculateZoneThreshold(self, dist_name = None, location = None, scale = None, zone_probability = None, 
                               sample = None, chromosome = None, param_type = None, 
                               reverse_transform = True):
        """
        Translate the zone probability to a threshold that can be compared against the signal 
        for which a distribution was fitted.

        params:
            dist_name:         Name of the fitted distribution, e.g. "Norm".
            location:          Position of distribution. Can be set instead of providing a param_type.
            scale:             Measure of distribution spread. Can be set instead of providing a 
                               param_type.
            zone_probability:  Probabilty threshold (between 0 to 1) from which the signal cut off is 
                               dervied from the distribution fitted for signal zone prediction.
            sample:            Sample name or sample ID to calculate zone threshold for.
            chromosome:        Chromosome that distribution was fitted to.
            param_type:        Type of parameter used when fitting the distribution, i.e. 'scipy_fit', 
                               'mean_fit' or 'median_fit'. Can be bypassed by setting custom location 
                               and scale.
            reverse_transform: Whether to reverse log transformation if applied when fitting 
                               distribution.
        """

        bw_idx = None

        if dist_name is None:
            dist_name = self.zone_distribution
        elif dist_name not in self.test_distributions:
            raise ValueError(f"Distribution {dist_name} was not found in test_distributions")
        
        if zone_probability is None:
            # Use the set threshold
            zone_probability = self.zone_probability
        elif (zone_probability <= 0) or (zone_probability >= 1):
            raise ValueError("zone_probability must be set as a probability between 0 and 1")

        if (location == None) and (scale == None):
            if (sample != None) and (chromosome != None):
                if isinstance(sample, str):
                    # Convert sample name to sample index
                    bw_idx = self.sampleToIndex(sample)
                elif isinstance(sample, (int, np.integer)):
                    # Sample index was given
                    bw_idx = sample
                else:
                    raise ValueError("Unsupported type for sample")

                if param_type is None:
                    param_type = self.zone_param_type
                
                # Get custom sample name
                sample_name = self.getSampleNames(return_custom = True)[bw_idx]

                # Open the fitted distribution statistics
                dist_stats_df = self.getDistributionStats()

                sample_found = sample_name in np.array(dist_stats_df["sample"])
                chrom_found = chromosome in np.array(dist_stats_df["chrom"], 
                                                     dtype = self.chromosomes.dtype)
                dist_found = dist_name in np.array(dist_stats_df["distribution"])
                params_found = param_type in np.array(dist_stats_df["param_type"])

                if not (sample_found & chrom_found & dist_found & params_found):
                    if self.verbose > 0:
                        print(f"Warning: {chromosome} {sample_name} with parameter type {param_type} "
                              f"was not found in the distribution statistics.")
                        print("Recreating combined distribution statistics.")

                    # Recreate the statistics
                    self.combineStats(stats_type = "distribution")
                    dist_stats_df = self.getDistributionStats()

                # Get fitted distribution for the sample and chromosome
                dist_stats_df = dist_stats_df[(dist_stats_df["sample"] == sample_name) &
                                              (dist_stats_df["chrom"] == chromosome) &
                                              (dist_stats_df["distribution"] == dist_name) &
                                              (dist_stats_df["param_type"] == param_type)]

                if len(dist_stats_df.index) == 0:
                    raise ValueError(f"{chromosome} {sample_name} with parameter "
                                     f"type {param_type} still missing.")

                # Set distribution parameters
                location = dist_stats_df["location"].iloc[0]
                scale = dist_stats_df["scale"].iloc[0]

            else:
                raise ValueError("If no location and scale are given, sample and chromosome "
                                 "must be given")

        elif location == None:
            raise ValueError("Location cannot be None is scale is set")
        elif scale == None:
            raise ValueError("Scale cannot be None is location is set")

        # Get the scipy distribution
        distribution = self.scipy_distributions[dist_name]
        # Calculate the threshold from the fitted distribution
        zone_threshold = distribution.ppf(zone_probability, loc = location, scale = scale)

        if self.log_transform and reverse_transform:
            # Reverse log transformation if applied
            zone_threshold *= np.e

        if sample is not None:
            if bw_idx is None:
                if isinstance(sample, str):
                    # Convert sample name to sample index
                    bw_idx = self.sampleToIndex(sample)
                elif isinstance(sample, (int, np.integer)):
                    # Sample index was given
                    bw_idx = sample
                else:
                    raise ValueError("Unsupported type for sample")

            if self.mean_genome_signals[bw_idx] < 0:
                # Flip sign if signal is negative
                zone_threshold *= -1
        elif np.all(self.mean_genome_signals < 0):
            # Assume signal is negative because all samples are negative
            zone_threshold *= -1

        elif np.any(self.mean_genome_signals < 0) and (self.verbose > 0):
            # Some samples are positive, and others are negative
            print("Warning: Sample was not given so cannot determine sign. Returning the absolute "
                  "value of the threshold")

        return zone_threshold

    def reverseZoneThreshold(self, dist_name, zone_threshold, location = None, scale = None, 
                             bw_idx = None, chromosome = None, param_type = None,):
        """
        Convert a zone threshold back to a probability.

        params:
            dist_name:         Name of the fitted distribution, e.g. "Norm".
            zone_threshold:    Threshold calculated from fitted distribution at a given zone probabilty.
            location:          Position of distribution. Can be set instead of providing a param_type.
            scale:             Measure of distribution spread. Can be set instead of providing a 
                               param_type.
            bw_idx:            Sample ID to calculate zone threshold for.
            chromosome:        Chromosome that distribution was fitted to.
            param_type:        Type of parameter used when fitting the distribution, i.e. 'scipy_fit', 
                               'mean_fit' or 'median_fit'. Can be bypassed by setting custom location 
                               and scale.
        """
        
        if dist_name not in self.test_distributions:
            raise ValueError(f"Distribution {dist_name} was not found in test_distributions")

        if (location == None) and (scale == None):
            if (bw_idx != None) and (chromosome != None):
                if param_type is None:
                    param_type = self.zone_param_type

                sample_name = self.getSampleNames(return_custom = True)[bw_idx]

                # Open the fitted distribution statistics
                dist_stats_df = self.getDistributionStats()

                sample_found = sample_name in np.array(dist_stats_df["sample"])
                chrom_found = chromosome in np.array(dist_stats_df["chrom"], 
                                                     dtype = self.chromosomes.dtype)
                dist_found = dist_name in np.array(dist_stats_df["distribution"])
                params_found = param_type in np.array(dist_stats_df["param_type"])

                if not (sample_found & chrom_found & dist_found & params_found):
                    if self.verbose > 0:
                        print(f"Warning: {chromosome} {sample_name} with parameter type {param_type} "
                              f"was not found in the distribution statistics.")
                        print("Recreating combined distribution statistics.")

                    # Recreate the statistics
                    self.combineStats(stats_type = "distribution")
                    dist_stats_df = self.getDistributionStats()

                # Get fitted distribution for the sample and chromosome
                dist_stats_df = dist_stats_df[(dist_stats_df["sample"] == sample_name) &
                                              (dist_stats_df["chrom"] == chromosome) &
                                              (dist_stats_df["distribution"] == dist_name) &
                                              (dist_stats_df["param_type"] == param_type)]

                if len(dist_stats_df.index) == 0:
                    raise ValueError(f"{chromosome} {sample_name} with parameter "
                                     f"type {param_type} still missing.")

                # Set distribution parameters
                location = dist_stats_df["location"].iloc[0]
                scale = dist_stats_df["scale"].iloc[0]

            else:
                raise ValueError("If no location and scale are given, bw_idx, chromosome "
                                 "must be given")

        elif location == None:
            raise ValueError("Location cannot be None is scale is set")
        elif scale == None:
            raise ValueError("Scale cannot be None is location is set")

        if self.log_transform:
            # Reverse log transformation if applied
            zone_threshold /= np.e

        # Get the scipy distribution
        distribution = self.scipy_distributions[dist_name]
        # Compute the probability
        zone_probability = distribution.cdf(zone_threshold, loc = location, scale = scale)

        return zone_probability

    def signalToZones(self, bw_idx, chromosome, dist_name, param_type, 
                      signal_type = "signal_transformed", zone_probability = None, extend_depth = 0, 
                      merge_depth = 0, min_region_bps = 35, quality_filter = True, 
                      min_different_bps = 5, chrom_fragment = 0, chrom_mean = 0, genome_mean = 0, 
                      round_to_bins = True):
        """
        Extract signal zone, i.e. coordinates of regions predicted to contain signal.

        params:
            bw_idx:                Sample ID to extract predicted regions of signal for.
            chromosome:            The chromosome to get signal zone coordinates for.
            dist_name:             Name of the distribution to fit from which zones are derived.
            param_type:            Type of parameters for fitting the distribution, i.e. 'scipy_fit', 
                                   'mean_fit' or 'median_fit'.
            signal_type:           Either 'transformed' to fit based on smoothed (and under 
                                   recommendation log transformed and downsampled) data, or 'signal' 
                                   to fit on the original signal.
            zone_probability:      Probabilty threshold (between 0 to 1) from which the signal cut off 
                                   is dervied from the distribution fitted for signal zone prediction.
            extend_depth:          Number of base pairs to add either side of each predicted signal zone.
            merge_depth:           Minimum distance between two zones for them to be combined into one.
            min_region_bps:        Minimum size of a region initially predicted to have signal, i.e. 
                                   the number of base pairs that must consecutively exceed the zone 
                                   threshold. The zone is likely to be larger than this if rounding 
                                   region coordinates to the nearest bins and adding one or more bins 
                                   worth of padding.
            quality_filter:        Set a True to filter out zones with insufficient signal.
            min_different_bps:     If using the quality filter, set the minimum number of different 
                                   base pairs that need to be found consecutively for a signal to be 
                                   considered good quality.
            chrom_fragment:        The closest value to zero of signal for the sample at the chromosome.
                                   Estimated magnitute of signal for a single fragment/read.
            chrom_mean:            The average signal (excluding zeros) for the sample's chromosome.
            genome_mean:           The average siganl (excluding zeros) across all the sample's 
                                   chromosomes.
            round_to_bins:         Whether to round region coordinates to the nearest bin coordinates,
                                   e.g. [74,756, 90,474] to [74,000, 91,000] for 1,000 base pair bins.
        """

        # Get the name of the sample according to its ID
        custom_sample_name = self.getSampleNames(return_custom = True)[bw_idx]
        file_sample_name = self.file_sample_names[bw_idx]
        chrom_size = self.chrom_sizes[chromosome]

        if param_type not in self.param_types:
            raise ValueError(f'Unknown parameter type "{param_type}". '
                             f"Options include {', '.join(self.param_types)}.")

        if zone_probability is None:
            zone_probability = self.zone_probability
        elif (zone_probability <= 0) or (zone_probability >= 1):
            raise ValueError("zone_probability must be set as a probability between 0 and 1")

        if self.verbose > 0:
            print(f"Predicting signal zones for {chromosome} {custom_sample_name} "
                  f"({bw_idx + 1}/{len(self.sample_names)})")

        if signal_type == "signal":
            # Open original signal
            signal = self.signalReader(bw_idx = bw_idx,
                                       start_idx = 0,
                                       end_idx = -1,
                                       chromosome = chromosome,
                                       pad_end = True,
                                       verbose = 0)
            
        elif signal_type == "signal_transformed":
            # Open smoothed signal
            chrom_signal_dir = os.path.join(self.output_directories["smooth_signal"], chromosome)
            signal = self.loadGzip(file_name = os.path.join(chrom_signal_dir, 
                                                            f"smooth-signal_{file_sample_name}.npy.gz"))

        else:
            raise ValueError(f"Unknown signal_type {signal_type}.\n"
                             f'signal_type must be set as either "transformed" to predict zones from '
                             f'the smoothed signal, or "signal" to predict from the original signal.')

        # Check if signal is positive or negative
        if self.mean_genome_signals[bw_idx] < 0:
            negative = True
        else:
            negative = False

        # Calculate the threshold from the fitted distribution
        zone_threshold = self.calculateZoneThreshold(sample = bw_idx,
                                                     chromosome = chromosome,
                                                     dist_name = dist_name, 
                                                     param_type = param_type,
                                                     zone_probability = zone_probability)

        # Get values which exceed the zone threshold (or reverse for negative)
        threshold_mask = np.abs(signal) > np.abs(zone_threshold)

        # Combine with mask that excludes positions with consecutive zeros
        signal_mask = self.createSignalMask(bw_idx = bw_idx,
                                            chromosome = chromosome,
                                            mask_size = len(signal))
        threshold_mask = threshold_mask & signal_mask

        # Exclude any masked base pair
        signal_idxs = np.where(threshold_mask)[0]
        del signal_mask
        del signal

        if (extend_depth > 0) or round_to_bins:
            pad_zones = True
        else:
            pad_zones = False

        if len(signal_idxs) == 0:
            if self.verbose > 0:
                print(f"No signal zones found for {chromosome} {custom_sample_name} as nothing "
                      f"exceeded the threshold of {np.round(zone_threshold, 3)}")

            # Set empty placeholders
            unpadded_signal_zones = np.empty((0,2), dtype = np.uint32)

            if pad_zones:
                padded_signal_zones = np.empty((0,2), dtype = np.uint32)

        else:
            # Find sequential indexes of signal and keep the start and end as coordinates
            # e.g. for [noise, signal, signal, signal, noise, noise] find [1,3]
            zone_coords = np.split(signal_idxs, np.where(np.diff(signal_idxs) != 1)[0] + 1)
            del signal_idxs
            zone_coords = np.array([np.array([pair[0], pair[-1]]) for pair in zone_coords], 
                                    dtype = np.uint32)
            # Filter to find those where the size exceeds a minimum number of base pairs
            zone_coords = zone_coords[np.where(np.diff(zone_coords) >= min_region_bps)[0]]

            if quality_filter:
                if negative:
                    # Flip stats and signal if signal is all negative
                    genome_mean *= -1
                    chrom_mean *= -1
                    chrom_fragment *= -1

                # Set a quality filter threshold from whichever is greater between:
                # 1) A scaled average signal across the genome
                # 2) A scaled average signal across the chromosome
                # 3) Signal equivalent to more than 10 overlapping fragments being detected
                min_threshold = max(genome_mean * 1.5,
                                    chrom_mean * 1.5,
                                    chrom_fragment * 10)

            else:
                # Allow all coordinates to pass quality control step if no filter applied
                passed_qc = True

            # Initialise list to store coordinates predicted to have signal prior to padding
            unpadded_signal_zones = []

            if pad_zones:
                # Record coordinates of signal after padding
                padded_signal_zones = []

            # Iterate over zones prior to filtering or padding
            for start, end in zone_coords:
                # Quality filtering only keep zones with a minimum number of consecutive base pairs 
                # over the quality threshold
                if quality_filter:
                    # Read signal within zone's coordinates
                    zone_signal = self.signalReader(bw_idx, chromosome, start_idx = start, end_idx = end)

                    if negative:
                        zone_signal *= -1

                    # Initialise counter
                    signal_count = 0
                    previous_value = 0
                    # Record whether the zone passed quality control
                    passed_qc = False

                    # Check values in the zone's signal until either part of the signal is determined as
                    # good quality or until the last value
                    for value in zone_signal:
                        if previous_value != value:
                            if value >= min_threshold:
                                # Increment the counter if a new value is found that is above the 
                                # noise threshold
                                signal_count += 1
                                if signal_count == min_different_bps:
                                    # If enough different values were found, the region passes quality 
                                    # control
                                    passed_qc = True
                                    break
                            else:
                                # Reset the counter as value was below the threshold
                                signal_count = 0
                            
                        # Keep record of the current value to check for consecutive duplication
                        previous_value = value

                if passed_qc:
                    # Convert types to prevent overflow warning
                    start = int(start)
                    end = int(end)

                    # Check if any unpadded zones recorded
                    if len(unpadded_signal_zones) == 0:
                        # Record the first
                        unpadded_signal_zones.append([start, end])
                    # Check if new start is within a set distance of the previous stop
                    elif start <= (unpadded_signal_zones[-1][1] + merge_depth):
                        # Found continuation so extend the previous zone's boundary
                        unpadded_signal_zones[-1][1] = end
                    # Check that the new zone is not equivalent to the previous one
                    elif end != unpadded_signal_zones[-1][1]:
                        # Record the new zone if different
                        unpadded_signal_zones.append([start, end])

                    if pad_zones:
                        if round_to_bins:
                            # Add padding and round start down to the nearest bin / round end up to 
                            # the nearest bin
                            start = int(max(0, np.floor((start - extend_depth) / 
                                                        self.bin_size)) * self.bin_size)
                            end = int(min(chrom_size, (np.ceil((end + extend_depth) / 
                                                               self.bin_size)) * self.bin_size))

                        else:
                            # Only add padding
                            start = max(0, start - extend_depth)
                            end = min(chrom_size, end + extend_depth)

                        # Record the zone after padding
                        if len(padded_signal_zones) == 0:
                            padded_signal_zones.append([start, end])
                        elif start <= (padded_signal_zones[-1][1] + merge_depth):
                            padded_signal_zones[-1][1] = end
                        elif end != padded_signal_zones[-1][1]:
                            padded_signal_zones.append([start, end])

            # Convert to numpy array
            unpadded_signal_zones = np.array(unpadded_signal_zones, dtype = np.uint32).reshape(-1, 2)

            if pad_zones:
                padded_signal_zones = np.array(padded_signal_zones, dtype = np.uint32).reshape(-1, 2)

            if self.verbose > 0:
                n_unpad_zones = unpadded_signal_zones.shape[0]
                if n_unpad_zones > 0:
                    n_pad_zones = padded_signal_zones.shape[0]
                    print(f"Keeping {n_unpad_zones} unpadded and {n_pad_zones} padded signal "
                          f"zone{'s' if n_unpad_zones != 1 else ''} for {chromosome} "
                          f"{custom_sample_name}")
                else:
                    print(f"No signal zones found for {chromosome} {custom_sample_name}")

        unpadded_zones_file = os.path.join(self.output_directories["signal_zones"], 
                                           chromosome, 
                                           "Unpadded",
                                           f"unpadded-zones_{file_sample_name}")

        # Save zones to file
        self.saveGzip(file_name = unpadded_zones_file,
                      array = unpadded_signal_zones)
        self.calculateZoneStats(bw_idx = bw_idx,
                                chromosome = chromosome,
                                use_padded = False,
                                use_merged = False, 
                                zones = unpadded_signal_zones,
                                exclude_zero = self.exclude_zero, 
                                zone_remove_percent = self.zone_remove_percent)
        
        if pad_zones:
            padded_zones_file = os.path.join(self.output_directories["signal_zones"],
                                             chromosome,
                                             "Padded",
                                             f"padded-zones_{file_sample_name}")

            self.saveGzip(file_name = padded_zones_file,
                          array = padded_signal_zones)
            self.calculateZoneStats(bw_idx = bw_idx,
                                    chromosome = chromosome,
                                    use_padded = True,
                                    use_merged = False, 
                                    zones = padded_signal_zones,
                                    exclude_zero = self.exclude_zero, 
                                    zone_remove_percent = self.zone_remove_percent)
            
    def mergeOverlapZones(self, chromosome = "", chrom_zones = np.array([]), sample_ids = np.array([]), 
                          padded_zones = False, merged_zones_file = "", verbose = 0):
        """
        Combine predicted regions of signal across samples.
        
        params:
            chromosome:         Chromosome to merge zones for (required if opening zones from file).
            chrom_zones:        Optionally provide array of pairs of start and end coordinates of zones 
                                for the chromosome. Otherwise zones are opened from file.
            sample_ids:         Array of sample indexes to merge regions for (use all samples by 
                                default).
            padded_zones:       Set as True to open padded zones, or False for unpadded zones.
            merged_zones_file:  If a file name is given, merged zones are saved to file. Otherwise the 
                                array of zones is returned.
            verbose:            Set at a value exceeding 0 to print progress/error messages.
        """

        if len(sample_ids) == 0:
            # Use all samples
            sample_ids = self.sample_ids

        # Check if coordinates provided
        if len(chrom_zones) > 0:
            custom_coords = True
        else:
            custom_coords = False

            if chromosome == "":
                raise ValueError(f"Chromosome must be specified to read "
                                 f"{'padded' if padded_zones else 'unpadded'} zones from file")

            if padded_zones:
                zones_dir = os.path.join(self.output_directories["signal_zones"], chromosome, "Padded")
                file_prefix = "padded-zones"
            else:
                zones_dir = os.path.join(self.output_directories["signal_zones"], chromosome, "Unpadded")
                file_prefix = "unpadded-zones"

            chrom_zones = []

            for bw_idx in sample_ids:
                # Open zones for the sample
                sample_zones = self.readArrayFile(bw_idx = bw_idx,
                                                  directory = zones_dir,
                                                  file_prefix = file_prefix)[1]
                
                if len(sample_zones) > 0:
                    chrom_zones.append(sample_zones)

        if len(chrom_zones) > 0:
            # Combine regions and sort by start position
            chrom_zones = np.unique(np.sort(np.vstack(chrom_zones), axis = 0), axis = 0)

            if verbose > 0:
                if custom_coords:
                    print("Merging coordinates")
                else:
                    print(f"Merging {'padded' if padded_zones else 'unpadded'} zones for {chromosome}")

            # Combine overlapping regions
            merged_zones = [chrom_zones[0]]

            for row_idx in range(1, len(chrom_zones)):
                # Check if new start is less than or equal to previous end
                if (chrom_zones[row_idx][0] <= merged_zones[-1][1]):
                    # Check if new end exceeds previous end
                    if chrom_zones[row_idx][1] > merged_zones[-1][1]:
                        # Extend the end coordinate of the current range
                        merged_zones[-1][1] = chrom_zones[row_idx][1]
                else:
                    # New range did not overlap so add it to the end
                    merged_zones.append(chrom_zones[row_idx])

            # Convert to array
            merged_zones = np.vstack(merged_zones)

            if verbose > 0:
                if custom_coords:
                    print(f"Merging created {len(merged_zones)} regions")
                else:
                    print(f"Merging created {len(merged_zones)} "
                          f"{'padded' if padded_zones else 'unpadded'} zones")

        else:
            # Set as empty array to record that no reliable signal was detected 
            # across samples for this chromosome
            merged_zones = np.empty((0,2), dtype = np.uint32)

            if verbose > 0:
                if custom_coords:
                    print("No merged coordinates found")
                else:
                    print(f"No merged {'padded' if padded_zones else 'unpadded'} "
                          f"zones found for {chromosome}... "
                          f"this could be due to low coverage across this chromosome.")
                    
        if merged_zones_file:
            # Save merged zones to file
            self.saveGzip(file_name = merged_zones_file,
                          array = merged_zones)
        else:
            # Return array
            return merged_zones

    def setDefaultBestParams(self, best_params):
        """
        If multiple parameters fitted equally well, select one.

        params:
            best_params: List of parameter types.
        """

        # Count the number of inferred best parameter types (usually only one)
        n_best_params = len(best_params)

        if n_best_params > 1:
            # Default to median as robust to outliers
            if "median_fit" in best_params:
                param_type = "median_fit"
            elif "scipy_fit" in best_params:
                param_type = "scipy_fit"
            else:
                # Set the type of parameters to use
                param_type = best_params[0]
        else:
            raise ValueError("best_params was empty")

        return param_type

    def predictSignalZones(self, chromosomes = [], zone_distribution = None, param_type = None,
                           extend_depth = None, merge_depth = None, round_to_bins = True, 
                           replace_existing = False):
        """
        Separates predicted regions of chromosome signal (signal zones) from noise per sample.

        params:
            chromosome:        List of chromosomes to predict signal zones for.
            zone_distribution: Set to force use of a specific distribution.
            param_type:        Set to force use of a specific parameter type for the distribution.
            extend_depth:      Number of base pairs to add either end of each predicted zone.
            merge_depth:       Minimum distance between two zones for them to be combined into one.
            round_to_bins:     Whether to round signal zone coordinates to the nearest full sized bins.
            replace_existing:  Whether to overwrite previously created files.
        """

        custom_sample_names = self.getSampleNames(return_custom = True)

        if len(chromosomes) == 0:
            chromosomes = self.chromosomes

        if len(self.sample_names) > 1:
            # Only create merged zones if processing more than one sample
            create_merged = True
        else:
            create_merged = False

        if extend_depth is None:
            extend_depth = self.extend_n_bins * self.bin_size
        if merge_depth is None:
            merge_depth = self.merge_depth
        else:
            self.setMergeDepth(merge_depth)

        if (extend_depth > 0) or round_to_bins:
            pad_zones = True
        else:
            pad_zones = False

        if zone_distribution is None:
            # Use preset distribution
            zone_distribution = self.zone_distribution

        if zone_distribution == "infer":
            if self.best_zone_distribution is None:
                # Predict the best distribution and parameter type
                best_zone_distribution = self.inferBestDistribution()
            else:
                best_zone_distribution = self.best_zone_distribution

            # Count the number of best distributions (one for the majority of cases)
            best_dist_names = list(best_zone_distribution.keys())
            n_best_dists = len(best_zone_distribution.keys())

            if n_best_dists > 1:
                print(f"{n_best_dists} best distributions were inferred.")
                print("Please manually select one by setting zone_distribution to continue.")
                print(f'Best distributions include: {", ".join(best_dist_names)}')
                return None
            else:
                # Set the name of the best distribution
                zone_distribution = list(best_zone_distribution.keys())[0]

            # Set the type of parameters to use
            param_type = self.setDefaultBestParams(best_zone_distribution[zone_distribution])
            
        elif param_type is None:
            # Use preset parameter type
            param_type = self.zone_param_type

        if param_type == "infer":
            # Using the given distribution, estimate the best parameter type
            best_zone_distribution = self.inferBestDistribution(zone_distribution = zone_distribution)
            param_type = self.setDefaultBestParams(best_zone_distribution[zone_distribution])

        # Create tuples of chromosomes and samples to process
        process_chrom_bws = []
        process_chroms = []

        for chrom in chromosomes:
            # Check if unpadded zones have been created
            chrom_unpadded_zones_dir = os.path.join(self.output_directories["signal_zones"],
                                                    chrom, "Unpadded")
            check_dirs = [chrom_unpadded_zones_dir]
            file_regex = "^unpadded-zones_*|"

            if pad_zones:
                # Also check for padded zones
                chrom_padded_zones_dir = os.path.join(self.output_directories["signal_zones"], 
                                                      chrom, "Padded")
                check_dirs.append(chrom_padded_zones_dir)
                file_regex += "^padded-zones_*|"

            file_regex += r"\.npy\.gz|"

            # Find indexes of samples where zones have not yet been saved to file
            incomplete_bw_idxs = self.findCurrentFiles(directories = check_dirs,
                                                       file_names_regex = file_regex,
                                                       replace_existing = replace_existing)
            
            if len(incomplete_bw_idxs) > 0:
                process_chroms.append(chrom)

                for bw_idx in incomplete_bw_idxs:
                    # Record that zones need to be prediced for the sample and chromosome
                    process_chrom_bws.append((chrom, bw_idx))

            for bw_idx in np.setdiff1d(self.sample_ids, incomplete_bw_idxs):
                # Also check that signal statistics exist for samples with zones
                stats_file = os.path.join(self.output_directories["signal_stats"], chrom, 
                                          f"signal-stats_{self.file_sample_names[bw_idx]}.csv")
                    
                if os.path.isfile(stats_file):
                    try:
                        # Attempt to read the statistics for the sample for a specific chromosome
                        stat_columns = pd.read_csv(stats_file)["signal_type"].tolist()

                        # Check if the expected statistics exist
                        missing_stats = (("signal_unpadded_sample_zones" not in stat_columns) or
                                         (pad_zones and ("signal_padded_sample_zones" not in stat_columns)) or
                                         (create_merged and ("signal_unpadded_merged_zones" not in stat_columns)) or
                                         (create_merged and pad_zones and ("signal_padded_merged_zones" not in stat_columns)))
                        
                        # Check if the expected statistics exist
                        if missing_stats:
                            process_chrom_bws.append((chrom, bw_idx))

                            if chrom not in process_chroms:
                                process_chroms.append(chrom)

                    except Exception as e:
                        print(f"Cannot find statistics for {custom_sample_names[bw_idx]}")
                        raise e
                else:
                    raise ValueError(f"Statistics file does not exist for {custom_sample_names[bw_idx]} "
                                     f"{chrom}")

        if self.interleave_sizes:
            process_chrom_bws = self.interleaveChroms(process_chrom_bws)

        if len(process_chrom_bws) > 0:
            if self.verbose > 0:
                print(f"Set to predict signal zones using {zone_distribution} distribution, "
                      f'{param_type.replace("_fit", "")} parameters and {round(self.zone_probability, 3)} '
                      f"zone probability")

            if self.quality_filter:
                chrom_fragments = {}
                chrom_means = {}
                global_means = {}

                for chrom in process_chroms:
                    # Get the statistics calculated for the unnormalised chromosome per sample
                    chrom_fragments[chrom], chrom_means[chrom], global_means[chrom] = self.getChromFragMean(chromosome = chrom, 
                                                                                                            sample_ids = self.sample_ids)

            if self.kernel is None:
                # Smoothing was not applied
                signal_type = "signal"
            else:
                # Predict based on smoothed signal
                signal_type = "signal_transformed"

            if (self.n_cores > 1) and (len(self.sample_ids) > 1):
                # Predict regions of signal per sample
                zone_processes = []

                with ProcessPoolExecutor(self.n_cores) as executor:
                    for chrom, bw_idx in process_chrom_bws:
                        zone_processes.append(executor.submit(self.signalToZones,
                                                              bw_idx = bw_idx,
                                                              chromosome = chrom,
                                                              dist_name = zone_distribution,
                                                              param_type = param_type,
                                                              signal_type = signal_type,
                                                              zone_probability = self.zone_probability,
                                                              extend_depth = extend_depth,
                                                              merge_depth = merge_depth,
                                                              min_region_bps = self.min_region_bps,
                                                              quality_filter = self.quality_filter,
                                                              min_different_bps = self.min_different_bps,
                                                              chrom_fragment = chrom_fragments[chrom][bw_idx] if self.quality_filter else 0,
                                                              chrom_mean = chrom_means[chrom][bw_idx] if self.quality_filter else 0,
                                                              genome_mean = global_means[chrom][bw_idx] if self.quality_filter else 0,
                                                              round_to_bins = round_to_bins))
                            
                    if self.checkParallelErrors(zone_processes):
                        return None

                    if create_merged:
                        # Release memory
                        del zone_processes
                        gc.collect()
                        merge_processes = []

                        # Merge zones across sample per chromosome
                        for chrom in process_chroms:
                            unpadded_merged_file = os.path.join(self.output_directories["signal_zones"], 
                                                                chrom, "Unpadded", "unpadded_merged_zones.npy.gz")
                            merge_processes.append(executor.submit(self.mergeOverlapZones,
                                                                   chromosome = chrom,
                                                                   padded_zones = False,
                                                                   merged_zones_file = unpadded_merged_file,
                                                                   verbose = self.verbose))

                            if pad_zones:
                                padded_merged_file = os.path.join(self.output_directories["signal_zones"], 
                                                                  chrom, "Padded", "padded_merged_zones.npy.gz")
                                merge_processes.append(executor.submit(self.mergeOverlapZones,
                                                                       chromosome = chrom,
                                                                       padded_zones = True,
                                                                       merged_zones_file = padded_merged_file,
                                                                       verbose = self.verbose))
                    
                        if self.checkParallelErrors(merge_processes):
                            return None

                        del merge_processes
                        gc.collect()
                        stats_processes = []

                        for chrom, bw_idx in process_chrom_bws:
                            # Open merged zone coordinates
                            merged_zones = self.getMergedZones(chrom, get_padded = False)["unpadded"]
                            
                            stats_processes.append(executor.submit(self.calculateZoneStats,
                                                                   bw_idx = bw_idx,
                                                                   chromosome = chrom,
                                                                   use_padded = False,
                                                                   use_merged = True,
                                                                   zones = merged_zones,
                                                                   exclude_zero = self.exclude_zero,
                                                                   zone_remove_percent = self.zone_remove_percent))

                        if self.checkParallelErrors(stats_processes):
                            return None

                        if pad_zones:
                            # Also save padded zone stats
                            gc.collect()
                            stats_processes = []

                            # Separate loop to prevent signal stats files being accessed simultaneously for the same sample
                            for chrom, bw_idx in process_chrom_bws:
                                merged_zones = self.getMergedZones(chrom, get_unpadded = False)["padded"]
                                stats_processes.append(executor.submit(self.calculateZoneStats,
                                                                       bw_idx = bw_idx,
                                                                       chromosome = chrom,
                                                                       use_padded = True,
                                                                       use_merged = True,
                                                                       zones = merged_zones,
                                                                       exclude_zero = self.exclude_zero,
                                                                       zone_remove_percent = self.zone_remove_percent))

                            if self.checkParallelErrors(stats_processes):
                                return None

            else:
                # Run sequentially rather than using parallelisation
                for chrom, bw_idx in process_chrom_bws:
                    self.signalToZones(bw_idx = bw_idx,
                                       chromosome = chrom,
                                       dist_name = zone_distribution,
                                       param_type = param_type,
                                       signal_type = signal_type,
                                       zone_probability = self.zone_probability,
                                       extend_depth = extend_depth,
                                       merge_depth = merge_depth,
                                       min_region_bps = self.min_region_bps,
                                       quality_filter = self.quality_filter,
                                       chrom_fragment = chrom_fragments[chrom][bw_idx] if self.quality_filter else 0,
                                       chrom_mean = chrom_means[chrom][bw_idx] if self.quality_filter else 0,
                                       genome_mean = global_means[chrom][bw_idx] if self.quality_filter else 0,
                                       round_to_bins = round_to_bins)

                if create_merged:
                    for chrom in process_chroms:
                        unpadded_merged_file = os.path.join(self.output_directories["signal_zones"], 
                                                            chrom,
                                                            "Unpadded",
                                                            "unpadded_merged_zones.npy.gz")
                        self.mergeOverlapZones(chromosome = chrom,
                                               padded_zones = False,
                                               merged_zones_file = unpadded_merged_file,
                                               verbose = self.verbose)

                        if pad_zones:
                            padded_merged_file = os.path.join(self.output_directories["signal_zones"], 
                                                              chrom,
                                                              "Padded",
                                                              "padded_merged_zones.npy.gz")
                            self.mergeOverlapZones(chromosome = chrom,
                                                   padded_zones = True,
                                                   merged_zones_file = padded_merged_file,
                                                   verbose = self.verbose)

                        for chrom, bw_idx in process_chrom_bws:
                            # Open merged zone coordinates
                            merged_zones = self.getMergedZones(chrom, get_padded = False)["unpadded"]
                            self.calculateZoneStats(bw_idx = bw_idx,
                                                    chromosome = chrom,
                                                    use_padded = False,
                                                    use_merged = True,
                                                    zones = merged_zones,
                                                    exclude_zero = self.exclude_zero,
                                                    zone_remove_percent = self.zone_remove_percent)

                            if pad_zones:
                                merged_zones = self.getMergedZones(chrom, get_unpadded = False)["padded"]
                                self.calculateZoneStats(bw_idx = bw_idx,
                                                        chromosome = chrom,
                                                        use_padded = True,
                                                        use_merged = True,
                                                        zones = merged_zones,
                                                        exclude_zero = self.exclude_zero,
                                                        zone_remove_percent = self.zone_remove_percent)

        elif self.verbose > 0:
            if create_merged:
                print(f"Signal and merged zone files already created for chromosomes")
            else:
                print(f"Signal zone files already created for chromosomes")

    def removePercentile(self, signal, remove_percent, exclude_zero = False, remove_upper = True, 
                         remove_lower = True, copy_signal = True):
        """
        Remove the lower and/or upper nth percentile of signal.

        params:
            signal:         Array of signal to process.
            remove_percent: Set as a value greater than zero to remove an upper and/or lower percentile 
                            from the signal. e.g. remove_percent = 100 removes the 100th percentiles.
            exclude_zero:   Set as True to remove zeros from the signal before percentile removal.
            remove_upper:   Set as True to remove the upper percentile.
            remove_lower:   Set as True to remove the lower percentile.
            copy_signal:    Whether to copy the numpy array. If False, the original array will be altered.
        """

        if (not remove_upper) and (not remove_lower):
            raise ValueError("remove_upper and remove_lower cannot both be set as False")

        if len(signal) == 0:
            return signal

        if copy_signal:
            # Copy to prevent changing the original signal
            signal = signal.copy()

        # Sort low to high by signal magnitude
        signal = np.sort(signal)
        # Index of lower and upper nth percentile
        percentile_pos = int(len(signal) / remove_percent)

        if not exclude_zero:
            # Count number of base pairs with signal greater than zero
            non_zero_count = np.count_nonzero(signal)

            # Check if removing the nth percentile of the signal will result in all zeros
            if non_zero_count <= percentile_pos:
                if self.verbose > 0:
                    print(f"Removing the {remove_percent}th "
                            "percentile results in all zeros")
                return None
                                                
        if remove_upper and remove_lower:
            # Remove lower and upper nth percentile
            signal = signal[percentile_pos:(len(signal) - percentile_pos)]
        elif remove_upper:
            # Only remove the upper nth percentile
            signal = signal[:(len(signal) - percentile_pos)]
        elif remove_lower:
            # Only remove the lower nth percentile
            signal = signal[percentile_pos:]

        return signal

    def calculateZoneStats(self, bw_idx, chromosome, use_padded, use_merged, exclude_zero, 
                           zone_remove_percent, zones):
        """
        Calculate signal statistics within zones for a specific sample and chromosome.
        
        params:
            bw_idx:              Sample ID to calculate statistics for.
            chromosome:          The chromosome to calculate statistics for.
            use_padded:          Set as True to calculate statistics for padded zones, or False for 
                                 unpadded zones.
            use_merged:          Set as True to calculate statistics for merged zones, or False for 
                                 sample specific zones.
            exclude_zero:        Set as True to remove zeros from signal when calculating statistics.
            zone_remove_percent: Set as a value greater than zero to remove an upper and/or lower 
                                 percentile from the signal. e.g. remove_percent = 100 removes the 
                                 100th percentiles.
            zones:               Array of pairs of zone start and end coordinates.
        """
            
        # Get the name of the sample according to its ID
        custom_sample_name = self.getSampleNames(return_custom = True)[bw_idx]
        file_sample_name = self.file_sample_names[bw_idx]
        n_samples = len(self.sample_ids)

        if self.verbose > 1:
            print(f"Calculating statistics across {'merged' if use_merged else 'sample'} "
                  f"{'padded' if use_padded else 'unpadded'} zones for {custom_sample_name} "
                  f"{chromosome} ({bw_idx + 1}/{n_samples})")

        # if len(zones) == 0:
        #     if use_merged:
        #         # Open merged zone coordinates
        #         if use_padded:
        #             zones = self.getMergedZones(chromosome, get_unpadded = False)["padded"]
        #         else:
        #             zones = self.getMergedZones(chromosome, get_padded = False)["unpadded"]
        #     else:
        #         # Open sample specific zone coordinates
        #         if use_padded:
        #             zones_dir = os.path.join(self.output_directories["signal_zones"], chromosome, "Padded")
        #             file_prefix = "padded-zones"
        #         else:
        #             zones_dir = os.path.join(self.output_directories["signal_zones"], chromosome, "Unpadded")
        #             file_prefix = "unpadded-zones"
                    
        #         zones = self.readArrayFile(bw_idx = bw_idx,
        #                                    directory = zones_dir,
        #                                    file_prefix = file_prefix)[1]

        # Count the number of base pairs across zones
        bp_coverage = int(np.sum(np.diff(zones)))

        if bp_coverage > 0:
            zones_found = True
            zone_signal = np.zeros(bp_coverage, dtype = np.float32)

            # Set first sample zone to check
            start_idx = 0
            end_idx = 0
        else:
            # No sample zone's to check
            zones_found = False
            zone_signal = np.array([], dtype = np.float32)

        if zones_found:
            for zone_idx in range(len(zones)):
                start_coord, end_coord = zones[zone_idx]
                zone_length = end_coord - start_coord
                end_idx = start_idx + zone_length
                zone_signal[start_idx:end_idx] = self.signalReader(bw_idx = bw_idx,
                                                                   start_idx = start_coord,
                                                                   end_idx = end_coord,
                                                                   chromosome = chromosome,
                                                                   verbose = 0)
                # Update for the next zone
                start_idx = end_idx

            if exclude_zero:
                # Remove zeros so they do not contribute to calculating statistics
                zone_signal = zone_signal[np.where(zone_signal != 0)]

            if zone_remove_percent > 0:
                # Remove the lower and upper nth percentile
                zone_signal = self.removePercentile(signal = zone_signal,
                                                    remove_percent = zone_remove_percent,
                                                    exclude_zero = exclude_zero,
                                                    copy_signal = False)

        # Calculate averages and deviations across the zones (or filler for empty signal)
        zone_stats = self.calculateAverages(signal = zone_signal)
        # Set file to save stats to
        signal_stats_file = os.path.join(self.output_directories["signal_stats"],
                                         chromosome, f"signal-stats_{file_sample_name}.csv")

        if self.verbose > 1:
            print(f"Saving statistics across zones for {custom_sample_name} {chromosome} "
                  f"({bw_idx + 1}/{n_samples})")

        signal_type = f"signal_{'padded' if use_padded else 'unpadded'}_"
        signal_type += f"{'merged' if use_merged else 'sample'}_zones"
        self.saveSignalStats(bw_idx = bw_idx,
                             chromosome = chromosome,
                             signal_stats = zone_stats,
                             signal_type = signal_type,
                             file_name = signal_stats_file)

    def normaliseChromSignal(self, bw_idx, chromosome, sample_norm_stats, file_name = ""):
        """
        Transform signal for a sample chromosome.
        
        params:
            bw_idx:            ID of the sample to normalise.
            chromosome:        Name of the chromosome, e.g. "chr1".
            sample_norm_stats: Dictionary of statistics to use for normalisation.
            file_name:         File name and path to save normalised signal to.
        """

        sample_name = np.array(self.getSampleNames(return_custom = True))[bw_idx]

        # Open full signal prior to normalisation
        signal = self.signalReader(bw_idx = bw_idx,
                                   start_idx = 0,
                                   end_idx = -1,
                                   chromosome = chromosome,
                                   verbose = 0)

        if self.verbose > 0:
            if self.norm_method == "Power":
                norm_message = f"Raising signal to the power of {self.norm_power} for"
            elif self.norm_method == "Log":
                norm_message = "Log transforming"
            else:
                norm_message = f"{self.norm_method} normalising"

            print(f"{norm_message} {chromosome} {sample_name} " 
                  f"({bw_idx + 1}/{len(self.sample_ids)})")
            
        if self.norm_method == "ZEN":
            sample_sd = sample_norm_stats["SD"]

            # Normalise the signal by variance
            signal = (signal / sample_sd).astype(np.float32)

        elif self.norm_method == "Power":
            # Raise signal to a power
            signal = (signal ** self.norm_power).astype(np.float32)

        elif self.norm_method == "Log":
            # Log transform signal
            signal = np.log(signal + 1).astype(np.float32)

        if file_name:
            if self.verbose > 0:
                if self.norm_method == "Power":
                    norm_message = f"signal raised to the power of {self.norm_power}"
                elif self.norm_method == "Log":
                    norm_message = "log transformed signal"
                else:
                    norm_message = f"{self.norm_method} normalised signal"

                print(f"Saving {norm_message} for {sample_name} {chromosome} " 
                      f"({bw_idx + 1}/{len(self.sample_ids)})")
                
            # Save to compressed numpy file
            self.saveGzip(file_name = file_name, array = signal)
        else:
            return signal
        
    def normaliseSignal(self, norm_stats_type = None, replace_existing = False):
        """
        Use an inbuilt normalisation method to rescale the bigWig signal.

        params:
            norm_stats_type:  Statistics type to use if using ZEN normalisation.
            replace_existing: Whether to overwrite previously created files.
        """

        if self.norm_method == "No normalisation":
            if self.verbose > 0:
                print("No normalisation method was set")
            # Do nothing
            return None

        elif not (self.norm_method in self.bigwig_norm_methods):
            if self.norm_method in self.bam_norm_methods:
                if self.verbose > 0:
                    print(f"{self.norm_method} normalisation is only supported to be run on BAM files")
                return None
                
            else:
                invalid_norm = self.norm_method
                # Disable normalisation if method not recognised
                self.norm_method = "No normalisation"
                # Warn user
                raise ValueError(f"Invalid bigWig normalisation method {invalid_norm}."
                                 f"\nDefaulting to no normalisation.")

        if norm_stats_type is None:
            norm_stats_type = self.norm_stats_type
        else:
            self.setNormStatsType(norm_stats_type)

        create_normalised_signal = True
        main_norm_dir = os.path.join(self.output_directories["norm_signal"], self.norm_method)
        norm_signal_bw_dirs = [os.path.join(main_norm_dir, chrom) for chrom in self.chromosomes]

        # Default to create signal for all samples
        incomplete_bw_idxs = self.sample_ids

        if not os.path.exists(main_norm_dir):
            # Create path if doesn't exist
            os.makedirs(main_norm_dir)

        elif not replace_existing:
            # Check if normalised signals already exist
            file_names_regex = "^" + self.norm_method + r"-norm-signal_*|\.npy\.gz"
            incomplete_bw_idxs = self.findCurrentFiles(directories = norm_signal_bw_dirs, 
                                                       file_names_regex = file_names_regex,
                                                       replace_existing = replace_existing)
                
            if len(incomplete_bw_idxs) == 0:
                create_normalised_signal = False
                if self.verbose > 0:
                    print(f"Full {self.norm_method} normalised signal already created")

        custom_sample_names = np.array(self.getSampleNames(return_custom = True))

        if create_normalised_signal:
            # Set-up for multiprocessing
            if self.n_cores > 1:
                process_parallel = True
                process_executor = ProcessPoolExecutor(self.n_cores)
            else:
                process_parallel = False

            if self.norm_method in ["Power", "Log"]:
                use_signal_stats = False

            else:
                use_signal_stats = True
                # Recreate signal statistics
                signal_stats = self.getSignalStats(replace_existing = True)

                # Get statistics for samples and chromosomes for which to create normalised signal
                signal_stats = signal_stats.loc[(np.isin(signal_stats["sample"], 
                                                         custom_sample_names[incomplete_bw_idxs])) & 
                                                (np.isin(signal_stats["chrom"], self.chromosomes)) &
                                                (signal_stats["signal_type"] == norm_stats_type) &
                                                (signal_stats["sum"] > 0)]

                # Check that statistics were calculated
                if len(signal_stats) == 0:
                    raise ValueError(f"Cannot perform {self.norm_method} normalisation as no "
                                     f"signal statistics were found for {norm_stats_type} for any "
                                     f"sample.\n"
                                     f"These can be created with extractZoneSignal.")

                # Check that all samples are present in signal statistics
                missing_samples = np.setdiff1d(np.unique(signal_stats["sample"]), 
                                                         custom_sample_names[incomplete_bw_idxs])

                if len(missing_samples) > 0:
                    raise ValueError(f'Missing rows for sample(s) {", ".join(missing_samples)} within '
                                     f"signal statistics.\n"
                                     f'Samples include: {", ".join(missing_samples)}.\n'
                                     f"Remaining statistics can be created with extractZoneSignal.")

                # Warn user if any chromosomes are absent
                missing_chroms = np.setdiff1d(np.unique(signal_stats["chrom"]), self.chromosomes)

                if len(missing_chroms) > 0:
                    if self.verbose > 0:
                        print(f"Warning: Missing rows for {len(missing_chroms)} chromosome(s) within "
                              f"signal statistics.\n"
                              f'Chromosomes include: {", ".join(missing_chroms)}.\n'
                              f"Remaining statistics can be created with extractZoneSignal.")

            if process_parallel:
                # Keep record of processes for error checking
                normalise_processes = []

            sample_norm_stats = {}

            for bw_idx in incomplete_bw_idxs:
                custom_sample_name = custom_sample_names[bw_idx]

                if use_signal_stats:
                    sample_stats = signal_stats.loc[(signal_stats["sample"] == custom_sample_name)]
                    sample_norm_stats = {}

                    if self.norm_method == "ZEN":
                        # Zone Z-Score
                        sample_norm_stats["SD"] = float(np.nanmedian(sample_stats["SD"]))

                for chrom_idx, chrom in enumerate(self.chromosomes):
                    norm_signal_bw_dir = norm_signal_bw_dirs[chrom_idx]
                    norm_file_name = os.path.join(norm_signal_bw_dir,
                                                  f"{self.norm_method}-norm-signal_{self.file_sample_names[bw_idx]}")

                    os.makedirs(norm_signal_bw_dir, exist_ok = True)

                    if process_parallel:
                        # Run each process to normalise each chromosome signal per sample
                        normalise_processes.append(process_executor.submit(self.normaliseChromSignal,
                                                                           bw_idx = bw_idx,
                                                                           chromosome = chrom,
                                                                           sample_norm_stats = sample_norm_stats,
                                                                           file_name = norm_file_name))
                    else:
                        # Normalise each chromosome signal per sample sequentially
                        self.normaliseChromSignal(bw_idx = bw_idx,
                                                  chromosome = chrom,
                                                  sample_norm_stats = sample_norm_stats,
                                                  file_name = norm_file_name)

            if process_parallel:
                # Check for errors during multiprocessing
                if self.checkParallelErrors(normalise_processes):
                    return None
                    
                # Release memory
                del normalise_processes
                gc.collect()

            create_bws = incomplete_bw_idxs

        else:
            create_bws = []

            for bw_idx in self.sample_ids:
                # Check if normalised bigWigs exist
                norm_bw_file = os.path.join(self.output_directories["bigwig"], self.norm_method,
                                            f"{custom_sample_names[bw_idx]}_{self.norm_method}.bw")

                if not os.path.isfile(norm_bw_file):
                    create_bws.append(bw_idx)

        if len(create_bws) > 0:
            # Save normalised signal to bigWig files
            self.createNormalisedBigWigs(sample_ids = create_bws)

    def extractBlacklistRegions(self, blacklist_file):
        """ 
        Record blacklist regions per chromosome.
        
        params:
            blacklist_file: File path to a zipped BED file of blacklist regions to exclude, 
                            e.g. "hg38-blacklist.v2.bed.gz".
        """

        blacklist_regions = {}

        if blacklist_file != "":
            # Open blacklist regions from BED file
            blacklist_df = pd.read_csv(blacklist_file, sep = "\t", header = None, usecols = range(3))

            # Convert regions to numpy arrays per chromosome
            for chrom in self.chromosomes:
                blacklist_regions[chrom] = np.array(blacklist_df.loc[blacklist_df[0] == chrom][[1,2]], 
                                                    dtype = np.uint32)

        return blacklist_regions

    def createMissingRegionFiles(self, blacklist_file = ""):
        """
        Save CSV files per sample containing potential deletions.
        
        params:
            blacklist_file: File path to a zipped BED file of blacklist regions to exclude, 
                            e.g. "hg38-blacklist.v2.bed.gz".
        """

        missing_regions_csv_dir = os.path.join(self.output_directories["missing_csv"], 
                                               "Missing_Regions")

        # Create empty directory to store intermediate files
        os.makedirs(missing_regions_csv_dir, exist_ok = True)

        if blacklist_file != "":
            if os.path.exists(blacklist_file):
                blacklist_regions = self.extractBlacklistRegions(blacklist_file)
            else:
                raise FileNotFoundError(f"Black list file {blacklist_file} does not exist")
        else:
            blacklist_regions = {}


        if self.n_cores > 1:
            with ProcessPoolExecutor(self.n_cores) as executor:
                # Run each process to read the chunked signal
                processes = [executor.submit(self.saveMissingRegions,
                                             bw_idx = bw_idx,
                                             output_dir = missing_regions_csv_dir,
                                             blacklist_regions = blacklist_regions) 
                                             for bw_idx in self.sample_ids]

                # Check for errors during multiprocessing
                if self.checkParallelErrors(processes):
                    return None
        else:
            # No parallelisation
            for bw_idx in self.sample_ids:
                self.saveMissingRegions(bw_idx = bw_idx, 
                                        output_dir = missing_regions_csv_dir,
                                        blacklist_regions = blacklist_regions)

    def saveMissingRegions(self, bw_idx, output_dir, blacklist_regions = {}):
        """
        Save all missing regions (potential deletions) for a sample to CSV.

        params:
            bw_idx:            Sample ID to save potential deletions for.
            blacklist_regions: Dictionary of chromosome coordinates to ignore.
            output_dir:        Directory to save CSV file to.
        """

        file_sample_name = self.file_sample_names[bw_idx]
        missing_dir = self.output_directories["missing_signal"]

        # Record all deletion coordinates across chromosomes in a single array
        missing_regions = []
        # Count number of rows corresponding to each chromosome
        chrom_rows = np.zeros((len(self.chromosomes)), dtype = np.uint32)

        for chrom_idx in range(len(self.chromosomes)):
            chrom = self.chromosomes[chrom_idx]

            # Open detected potential deletions
            regions_file = os.path.join(missing_dir, chrom,
                                        f"missing-regions_{file_sample_name}.npy.gz")
            chrom_missing_regions = self.loadGzip(file_name = regions_file)
                
            if chrom in blacklist_regions:
                # Check each blacklisted region for the chromosome
                chrom_blacklist_regions = blacklist_regions[chrom]
                region_mask = np.ones(len(chrom_missing_regions), dtype = bool)

                for region in chrom_blacklist_regions:
                    # Mask any regions that overlap a blacklisted region
                    region_mask &= ~((chrom_missing_regions[:,0] <= region[0]) & 
                                     (chrom_missing_regions[:,1] >= region[1]))

                chrom_missing_regions = chrom_missing_regions[region_mask]

            missing_regions.append(chrom_missing_regions)
            chrom_rows[chrom_idx] = len(chrom_missing_regions)

        missing_regions = np.vstack(missing_regions, dtype = np.uint32)

        chrom_col, start_col, end_col = self.region_coords_cols
        missing_df = pd.DataFrame({chrom_col: np.repeat(self.chromosomes, repeats = chrom_rows),
                                   start_col: missing_regions[:,0],
                                   end_col: missing_regions[:,1]})

        missing_df.to_csv(os.path.join(output_dir, 
                                       f"missing-regions_{self.getSampleNames(return_custom = True)[bw_idx]}.csv"), 
                          header = True, index = False)

    def saveSignalStats(self, bw_idx, chromosome, signal_stats, signal_type, file_name = "", 
                        return_df = False):
        """
        Save signal statistics to CSV for a sample.

        params:
            bw_idx:       Sample ID statistics were calculated for.
            chromosome:   Chromosome statistics were calculated for.
            signal_stats: Dictionary of statistics for the signal to save.
            signal_type:  Name for the new rows of signal statistics.
            file_name:    Optionally set this as a custom file name to save statistics to.
            return_df:    Set as True to return the statistics DataFrame.
        """

        # Get the name of the sample according to its ID
        custom_sample_name = self.getSampleNames(return_custom = True)[bw_idx]
        file_sample_name = self.file_sample_names[bw_idx]

        if file_name == "":
            file_name = os.path.join(self.output_directories["signal_stats"],
                                     chromosome,
                                     f"signal-stats_{file_sample_name}.csv")
        elif not file_name.endswith(".csv"):
            # Add numpy and gzip file extension
            file_name += ".csv"

        if signal_type == "signal":
            # Convert signal statistics dictionary to a dataframe
            stats_df = pd.DataFrame({"signal_type": [signal_type],
                                     "fragmentEstimate": [signal_stats["fragmentEstimate"]],
                                     "coverage": [signal_stats["coverage"]],
                                     "sum": [signal_stats["sum"]],
                                     "mean": [signal_stats["mean"]],
                                     "median": [signal_stats["median"]],
                                     "SD": [signal_stats["SD"]],
                                     "MAD": [signal_stats["MAD"]],
                                     "meanAD": [signal_stats["meanAD"]]})
            # Cast the number of base pairs covered as an integer
            stats_df["coverage"] = stats_df["coverage"].astype(int)
                                    
        else:
            # Check that the file exists
            if not os.path.isfile(file_name):
                raise ValueError(f"Signal statistics file does not exist for {custom_sample_name} "
                                 f"{chromosome}")
                
            # Check if the file is empty
            if os.path.getsize(file_name) == 0:
                # Delete an empty file as creation was likely interrupted
                os.remove(file_name)

                # Read original signal
                signal = self.signalReader(bw_idx = bw_idx,
                                           start_idx = 0,
                                           end_idx = -1,
                                           chromosome = chromosome,
                                           pad_end = True,
                                           verbose = 0)

                # Attempt to recreate the signal statistics for the original signal
                signal_stats = self.calculateAverages(signal = signal)
                self.saveSignalStats(bw_idx, chromosome, signal_stats, signal_type, file_name)

            if os.path.getsize(file_name) == 0:
                raise ValueError(f"Signal statistics for {custom_sample_name} {chromosome} are "
                                 f"missing.\n")

            try:
                # Read non-empty file
                stats_df = pd.read_csv(file_name)
            except Exception as e:
                # Catch unexpected errors
                print(f"Could not read signal statistics for {custom_sample_name} {chromosome} "
                      f'from file "{file_name}" due to exception:')
                raise e

            # Check if entry already exists for the signal type
            if signal_type in np.array(stats_df["signal_type"]):
                if self.verbose > 1:
                    print(f"Warning: overwriting {signal_type} row for {custom_sample_name} "
                          f"{chromosome}")

                row_idx = stats_df.loc[stats_df["signal_type"] == signal_type].index[0]
                
                # Add values in the signal stats dictionary
                for col in signal_stats.keys():
                    stats_df.loc[row_idx, col] = signal_stats[col]

            else:
                # Create an empty row to add to the dataframe
                new_row = pd.DataFrame([np.nan] * len(stats_df.columns)).T
                new_row.columns = stats_df.columns
                new_row["signal_type"] = signal_type

                # Add values in the signal stats dictionary
                for col in signal_stats.keys():
                    new_row[col] = signal_stats[col]

                # Add the new row
                stats_df = pd.concat([stats_df, new_row], ignore_index = True)

        # Write signal stats to CSV
        stats_df.to_csv(file_name, header = True, index = False)

        if return_df:
            return stats_df

    def saveDistributionStats(self, distribution_stats, file_name = ""):
        """
        Save statistics from distribution(s) fitted to smoothed signal to CSV.

        params:
            distribution_stats: Dictionary of statistic names and values.
            file_name:          Set as a file name to save to, else statistics are 
                                returned as a DataFrame.
        """

        # Get names of successfully fitted distributions
        dist_names = set(distribution_stats.keys())
        dist_names.discard("signal")
        dist_names.discard("signal_transformed")
        dist_names = np.array(list(dist_names))

        # Count the number of rows and identify columns
        n_rows = 0
        columns = np.array([])

        for dist_name in dist_names:
            dist_stats = distribution_stats[dist_name]
            n_rows += len(dist_stats.keys())
            for param_type in distribution_stats[dist_name]:
                param_stats = dist_stats[param_type]
                param_keys = np.array(list(param_stats.keys()))
                columns = np.unique(np.concatenate((columns, param_keys)))

        # Add the distribution name and parameter type columns
        columns = np.concatenate((["distribution", "param_type"], columns))

        # Initialise the empty dataframe of NaNs for distribution stats
        dist_df = pd.DataFrame(index = np.arange(n_rows), columns = columns, dtype = float)
        # Cast first two columns as strings
        dist_df["distribution"] = dist_df["distribution"].astype(str)
        dist_df["param_type"] = dist_df["param_type"].astype(str)

        row_idx = 0

        for dist_idx, dist_name in enumerate(dist_names):
            dist_stats = distribution_stats[dist_name]
            n_param_types = len(dist_stats.keys())
            row_idxs = np.arange(dist_idx * n_param_types, (dist_idx + 1)* n_param_types)
            
            # Set the distribution name and parameter type
            dist_df.iloc[row_idxs, 0] = dist_name

            for param_idx, param_type in enumerate(dist_stats):
                param_stats = dist_stats[param_type]
                dist_df.iloc[row_idx, 1] = param_type

                # Set the values
                for col, value in zip(param_stats.keys(), param_stats.values()):
                    dist_df.at[row_idx, col] = value

                row_idx += 1

        if file_name:
            if not file_name.endswith(".csv"):
                # Add numpy and gzip file extension
                file_name += ".csv"

            # Write distribution stats to CSV
            dist_df.to_csv(file_name, header = True, index = False)
        else:
            return dist_df

    def combineStats(self, stats_type, file_name = "", lock = None, return_df = False):
        """
        Combine statistics across samples into a CSV.

        params:
            stats_type: Set as either "signal" or "distribution".
            file_name:  Can specify a file name to read statistics CSV from. Otherwise uses default file.
            lock:       Can be set as a threading/multiprocessing lock to prevent parallel access to 
                        the statistics.
            return_df:  Set as True to return the updated statistics DataFrame.
        """

        custom_sample_names = np.array(self.getSampleNames(return_custom = True))
        file_sample_names = self.file_sample_names

        if stats_type not in ["signal", "distribution"]:
            raise ValueError('stats_type must be set as either "signal" or "distribution"')

        if file_name == "" and (not return_df):
            if stats_type == "signal":
                # Set directory and part of name of file to output
                output_directory = self.output_directories["output_stats"]
                output_file = "combined-signal-stats.csv"
            else:
                output_directory = self.output_directories["output_stats"]
                output_file = "combined-distribution-stats.csv"

            file_name = os.path.join(output_directory, output_file)
        else:
            output_directory = os.path.dirname(file_name)

        if output_directory:
            # Create directory if does not yet exist
            os.makedirs(output_directory, exist_ok = True)

        if stats_type == "signal":
            input_directory = self.output_directories["signal_stats"]
            file_prefix = "signal-stats"
        else:
            input_directory = self.output_directories["dist_stats"]
            file_prefix = "distribution-stats"

        # Record the statistics per sample per chromosome
        combined_stats = []

        for bw_idx in self.sample_ids:
            file_sample_name = file_sample_names[bw_idx]
            custom_sample_name = custom_sample_names[bw_idx]

            for chrom in self.chromosomes:
                stats_file = os.path.join(input_directory, chrom, f"{file_prefix}_{file_sample_name}.csv")
                
                if os.path.isfile(stats_file):
                    try:
                        # Attempt to read the statistics for the sample for a specific chromosome
                        stats_df = pd.read_csv(stats_file)
                        n_rows = len(stats_df)
                        # Add columns for sample name and chromosome
                        stats_df.insert(0, "chrom", [chrom] * n_rows)
                        stats_df.insert(0, "sample", [custom_sample_name] * n_rows)
                        combined_stats.append(stats_df)

                    except Exception as e:
                        print(f"Could not read {stats_type} statistics for {custom_sample_name} "
                              f"{chrom} due to exception:")
                        raise e
                else:
                    if self.verbose > 0:
                        print(f"Warning: Could not read {stats_type} statistics for "
                              f"{custom_sample_name} {chrom} as file does not exist")

        if len(combined_stats) == 0:
            if stats_type == "signal":
                raise ValueError(f"Could not find signal statistics for any chromosome.\n"
                                 f"These can be created with convolveSignals.")
            else:
                raise ValueError(f"Could not find distribution statistics for any chromosome.\n"
                                 f"To produce distribution statistics, run testDistributions.")

        # Combine list of dataframes
        combined_stats_df = pd.concat(combined_stats)

        if file_name:
            if lock is not None:
                with lock:
                    combined_stats_df.to_csv(file_name, header = True, index = False)
            else:
                combined_stats_df.to_csv(file_name, header = True, index = False)

        if return_df:
            return combined_stats_df

    def combineDistributionStats(self, file_name = "", return_df = False):
        """
        Combine distribution statistics across samples into a CSV.

        params:
            file_name: Can specify a file name to read statistics CSV from. Otherwise uses default file.
            return_df: Set as True to return the updated statistics DataFrame.
        """

        custom_sample_names = np.array(self.getSampleNames(return_custom = True))
        file_sample_names = self.file_sample_names

        if file_name == "" and (not return_df):
            file_name = os.path.join(self.output_directories["output_stats"], 
                                    "combined-distribution-stats.csv")
            # Make directory if does not exist
            os.makedirs(self.output_directories["output_stats"], exist_ok = True)

        # Record the distribution statistics per sample per chromosome
        combined_dist_stats = []

        for bw_idx in self.sample_ids:
            file_sample_name = file_sample_names[bw_idx]
            custom_sample_name = custom_sample_names[bw_idx]

            for chrom in self.chromosomes:
                dist_stats_file = os.path.join(self.output_directories["dist_stats"], chrom, 
                                               f"distribution-stats_{file_sample_name}.csv")
                
                if os.path.isfile(dist_stats_file):
                    try:
                        # Attempt to read the distribution statistics for the sample for a specific 
                        # chromosome
                        dist_df = pd.read_csv(dist_stats_file)
                        n_rows = len(dist_df)
                        # Add columns for sample name and chromosome
                        dist_df.insert(0, "chrom", [chrom] * n_rows)
                        dist_df.insert(0, "sample", [custom_sample_name] * n_rows)
                        combined_dist_stats.append(dist_df)

                    except Exception as e:
                        print(f"Could not read distribution statistics for {custom_sample_name} "
                              f"{chrom} due to exception:")
                        raise e
                else:
                    if self.verbose > 0:
                        print(f"Warning: Could not read distribution statistics for "
                              f"{custom_sample_name} {chrom} as file does not exist")

        # Combine list of dataframes
        combined_dist_stats_df = pd.concat(combined_dist_stats)

        if file_name:
            combined_dist_stats_df.to_csv(file_name, header = True, index = False)

        if return_df:
            return combined_dist_stats_df

    def createSmoothBigWigs(self, sample_ids = []):
        """
        Save signal smoothed via convolution to bigWig files.

        params:
            sample_ids: Same IDs to create bigWigs for.
        """

        if len(sample_ids) == 0:
            sample_ids = self.sample_ids
        elif max(sample_ids) > len(self.sample_ids):
            raise ValueError("Invalid sample IDs")

        smooth_bigwigs_dir = self.output_directories["smooth_bigwig"]
            
        # Create new directory to store smoothed bigwig tracks
        os.makedirs(smooth_bigwigs_dir, exist_ok = True)

        sample_names = np.array(self.getSampleNames(return_custom = True))[sample_ids]
                
        if self.n_cores > 1:
            # Run multiple processes in parallel
            processes = []
                
            with ProcessPoolExecutor(self.n_cores) as executor:
                for bw_idx, name in zip(self.sample_ids, sample_names):
                    processes.append(executor.submit(self.saveBigWig,
                                                     bw_idx = bw_idx,
                                                     file_name = f"{name}_smooth.bw",
                                                     directory = smooth_bigwigs_dir,
                                                     signal_type = "smooth"))
                self.checkParallelErrors(processes)

        else:
            for bw_idx, name in zip(self.sample_ids, sample_names):
                self.saveBigWig(bw_idx = bw_idx,
                                file_name = f"{name}_smooth.bw",
                                directory = smooth_bigwigs_dir,
                                signal_type = "smooth")

    def createNormalisedBigWigs(self, sample_ids = []):
        """
        Save normalised signals to bigWig files.

        params:
            sample_ids: Same IDs to create bigWigs for.
        """

        if len(sample_ids) == 0:
            sample_ids = self.sample_ids
        elif max(sample_ids) > len(self.sample_ids):
            raise ValueError("Invalid sample IDs")

        if self.norm_method != "No normalisation":
            normalised_bigwigs_dir = os.path.join(self.output_directories["bigwig"], self.norm_method)
            
            # Create new directory to store normalised bigwig tracks
            os.makedirs(normalised_bigwigs_dir, exist_ok = True)

            sample_names = np.array(self.getSampleNames(return_custom = True))[sample_ids]
                
            if self.n_cores > 1:
                # Run multiple processes in parallel
                processes = []
                
                with ProcessPoolExecutor(self.n_cores) as executor:
                    for bw_idx, name in zip(self.sample_ids, sample_names):
                        processes.append(executor.submit(self.saveBigWig,
                                                         bw_idx = bw_idx,
                                                         file_name = f"{name}_{self.norm_method}.bw",
                                                         directory = normalised_bigwigs_dir,
                                                         signal_type = "norm"))
                    self.checkParallelErrors(processes)

            else:
                for bw_idx, name in zip(self.sample_ids, sample_names):
                    self.saveBigWig(bw_idx = bw_idx,
                                    file_name = f"{name}_{self.norm_method}.bw",
                                    directory = normalised_bigwigs_dir,
                                    signal_type = "norm")
        else:
            print(f"Cannot create new normalised bigWigs for {len(sample_ids)} "
                  f"samples as no normalisation was applied")

    def createDistributionPlotValues(self, bw_idx, chromosome, signal = np.array([]), 
                                     transform_signal = True, downsample = False, 
                                     plot_distributions = [], calculate_stats = False):
        """
        Create a dictionary of distributions, signal and statistics for plotting fitted distributions.
        
        params:
            bw_idx:             The sample ID to get plot values for.
            chromosome:         The chromosome to get plot values for.
            signal:             Optionally provide a signal to optionally transform and then create 
                                the plot across.
            transform_signal:   Whether to transform the signal, e.g. by log transformation and/or 
                                downsampling.
            downsample:         Set as True to downsample the signal.
            plot_distributions: List of distributions to plot.
            calculate_stats:    Whether to include dictionary of signal statistics in output.
        """

        sample_name = self.getSampleNames(return_custom = True)[bw_idx]

        # Open the fitted distribution statistics
        dist_stats_file = os.path.join(self.output_directories["dist_stats"],
                                       chromosome,
                                       f"distribution-stats_{self.file_sample_names[bw_idx]}.csv")

        if not os.path.isfile(dist_stats_file):
            raise FileNotFoundError(f"Distribution statistics not found for {chromosome} {sample_name}")

        dist_stats_df = pd.read_csv(dist_stats_file)
        all_distributions = np.unique(dist_stats_df["distribution"])

        if len(plot_distributions) == 0:
            # If none specified, set to plot all distributions
            plot_distributions = all_distributions
        elif not (isinstance(plot_distributions, (list, np.ndarray))):
            raise ValueError("plot_distributions must be set as a list or array")
        else:
            keep_distributions = []

            for dist_name in plot_distributions:
                dist_name = dist_name.lower()
                if dist_name.startswith("norm"):
                    dist_name = "norm"
                keep_distributions.append(dist_name)

            plot_distributions = np.array(keep_distributions)
            missing_distributions = np.setdiff1d(plot_distributions, all_distributions)

            if len(missing_distributions) > 0:
                raise ValueError(f"Cannot plot distributions the following distributions as no "
                                f"statistics were calculated for them: "
                                f'"{", ".join(missing_distributions)}"')

            if len(plot_distributions) < len(all_distributions):
                # Subset specific distributions
                dist_stats_df = dist_stats_df[dist_stats_df["distribution"].isin(plot_distributions)]

        if len(signal) == 0:
            if self.kernel is None:
                signal = self.signalReader(bw_idx = bw_idx,
                                           start_idx = 0,
                                           end_idx = -1,
                                           chromosome = chromosome,
                                           pad_end = True,
                                           verbose = 0)
            else:
                # Read signal from array if not directly given
                chrom_signal_dir = os.path.join(self.output_directories["smooth_signal"], chromosome)
                signal = self.readArrayFile(bw_idx = bw_idx, directory = chrom_signal_dir, 
                                            file_prefix = "smooth-signal")[1]

        if transform_signal:
            if downsample:
                # Plot downsampled signal in the histogram
                downsample_size = self.downsample_size
            else:
                # Plot all signal in the histogram
                downsample_size = 0

            if calculate_stats:
                signal_to_fit, signal_stats = self.transformSignal(signal,
                                                                   log_transform = self.log_transform,
                                                                   downsample_size = downsample_size,
                                                                   calculate_stats = True)
            else:
                signal_to_fit = self.transformSignal(signal,
                                                     log_transform = self.log_transform,
                                                     downsample_size = downsample_size,
                                                     calculate_stats = False)
            del signal
        else:
            signal_to_fit = signal

            if calculate_stats:
                # Calculate statistics to plot for the signal
                mean = np.mean(signal_to_fit)
                median = np.median(signal_to_fit)
                sd = np.std(signal_to_fit)
                mad = np.median(np.abs(signal_to_fit - median))

                signal_stats = {"mean": mean,
                                "median": median,
                                "SD": sd,
                                "MAD": mad}

        results = {"plot_distributions": plot_distributions,
                   "signal_to_fit": signal_to_fit,
                   "dist_stats_df": dist_stats_df}
        
        if calculate_stats:
            results["signal_stats"] = signal_stats

        return results

    def plotDistributionFit(self, plot_sample, chromosome, signal = np.array([]), transform_signal = True,
                            downsample = False, plot_distributions = [], plot_param_types = [], 
                            n_bins = 30, plot_zone_threshold = True, zone_probability = None, 
                            title = "", max_ncols = 3, x_limits = [], y_limits = [],
                            plot_width = 16, plot_height = 4, pdf_name = ""):
        """
        Create histogram plot(s) of distribution(s) fitted for signal zone prediction. 
        These can help visually assess goodness-of-fit.

        params:
            plot_sample:         Name of sample to plot values for.
            chromosome:          The chromosome to get plot values for.
            signal:              Optionally provide a signal to optionally transform and then create 
                                 the plot across.
            transform_signal:    Whether to transform the signal, e.g. by log transformation and/or 
                                 downsampling.
            downsample:          Set as True to downsample the signal.
            plot_distributions:  List of distributions to plot.
            plot_param_types:    List of distribution fitting parameter types to plot, 
                                 e.g. ["mean_fit", "median_fit", "scipy_fit"].
            n_bins:              Number of bins for the histogram.
            plot_zone_threshold: Set as True to include a line per histogram at the value of the zone 
                                 threshold.
            zone_probability:    Probabilty threshold (between 0 to 1) from which the signal cut off 
                                 is dervied from the distribution fitted for signal zone prediction.
            title:               Title to display at top of plot.
            max_ncols:           Maximum number of columns for the sub-plots. If number of sub-plots 
                                 exceeds this, sub-plots are shown on separate rows.
            x_limits:            If set with two x-coordinates, e.g. [-10, 10], these will be used as 
                                 the x-axis range per histogram.
            y_limits:            If set with two y-coordinates, e.g. [0, 10], these will be used as 
                                 the x-axis range per histogram.
            plot_width:          Length of the plot.
            plot_height:         Height of the plot.
            pdf_name:            To save plot to PDF, set this as a file name.
        """

        custom_sample_names = self.getSampleNames(return_custom = True)

        if isinstance(plot_sample, str):
            bw_idx = self.sampleToIndex(plot_sample)

        elif isinstance(plot_sample, (int, np.integer)):
            # Use as sample ID
            bw_idx = plot_sample

            if bw_idx < 0 or bw_idx >= len(custom_sample_names):
                raise IndexError(f"Plot sample was set as a sample ID ({bw_idx}) that is out of range")
            
            plot_sample = custom_sample_names[bw_idx]

        sample_name = np.array(self.getSampleNames(return_custom = True))[bw_idx]
        plot_param_types = np.array(plot_param_types)

        if pdf_name:
            if not pdf_name.endswith(".pdf"):
                pdf_name = pdf_name + ".pdf"

            # Create directory to store plots
            os.makedirs(self.output_directories["plots"], exist_ok = True)

        plot_values = self.createDistributionPlotValues(bw_idx = bw_idx,
                                                        chromosome = chromosome,
                                                        signal = signal,
                                                        transform_signal = transform_signal,
                                                        downsample = downsample,
                                                        plot_distributions = plot_distributions,
                                                        calculate_stats = True)

        # Sanitised input distributions
        plot_distributions = plot_values["plot_distributions"]
        # The signal to plot a histogram of
        signal_to_fit = plot_values["signal_to_fit"]
        # Distribution statistics for the sample's chromosome
        dist_stats_df = plot_values["dist_stats_df"]
        # Contains mean, median, SD and MAD
        signal_stats = plot_values["signal_stats"]
        del plot_values

        # Find the unique parameter combinations
        param_types = np.sort(np.unique(dist_stats_df["param_type"]))

        if len(plot_param_types) > 0:
            # Find specific parameter types
            param_types = np.intersect1d(param_types, plot_param_types)

            if len(param_types) == 0:
                raise ValueError("No rows in distribution statistics were found for any value in "
                                 "plot_param_types")

        # Create dictionary to map each parameter type to its own axis
        n_param_types = len(param_types)
        param_types_ax_idxs = dict(zip(param_types, np.arange(n_param_types)))

        max_ncols = max(1, min(max_ncols, n_param_types))
        nrows = n_param_types // max_ncols

        if n_param_types % max_ncols != 0:
            nrows += 1

        # Create the subpplots for fitted distributions
        fig, axes = plt.subplots(nrows, max_ncols, figsize = (plot_width, plot_height))
        x = np.linspace(min(signal_to_fit), max(signal_to_fit), 1000)

        # Set colours for plotting distributions
        dist_colours = {"norm": "#1f77b4",
                        "laplace": "#ff7f0e",
                        "logistic": "#2ca02c",
                        "gumbel_l": "#9467bd",
                        "gumbel_r": "#e377c2"}

        if plot_zone_threshold:
            if zone_probability is None:
                zone_probability = self.zone_probability

        if len(x_limits) == 2:
            set_x_lims = True
        else:
            set_x_lims = False
        if len(y_limits) == 2:
            set_y_lims = True
        else:
            set_y_lims = False

        # Set colours for statistic lines
        average_colour = "black"
        deviation_colour = "red"
        # Line transparency
        line_alpha = 0.8

        alt_dist_titles = {"norm": "Normal",
                           "gumbel_l": "Gumbel Left",
                           "gumbel_r": "Gumbel Right"}

        axes = np.atleast_1d(axes).flatten()

        for ax_idx, ax in enumerate(axes):
            # Add the signal density
            ax.hist(signal_to_fit, bins = n_bins, density = True, color = "whitesmoke", 
                    edgecolor = "black", alpha = 0.4)

            # Set the title as the parameter type
            param_type = param_types[ax_idx]
            ax.set_title(param_type.replace("_", " ").title())

            if set_x_lims:
                ax.set_xlim(x_limits[0], x_limits[1])
            if set_y_lims:
                ax.set_ylim(y_limits[0], y_limits[1])

            if ax_idx == 0:
                # Set shared y-axis label
                ax.set_ylabel("Density")
            else:
                # Hide the y-axis label
                ax.tick_params(left = False, labelleft = False)

            if param_type == "mean_fit":
                ax.axvline(x = signal_stats["mean"], color = average_colour, 
                           linestyle = "-", alpha = line_alpha,
                           label = f'Mean = {round(signal_stats["mean"], 3)}')
                ax.axvline(x = signal_stats["SD"], color = deviation_colour, 
                           linestyle = "-", alpha = line_alpha,
                           label = f'SD = {round(signal_stats["SD"], 3)}')
                ax.legend()
            elif param_type == "median_fit":
                ax.axvline(x = signal_stats["median"], color = average_colour, 
                           alpha = line_alpha, linestyle = "-", 
                           label = f'Median = {round(signal_stats["median"], 3)}')
                ax.axvline(x = signal_stats["MAD"], color = deviation_colour, 
                           alpha = line_alpha, linestyle = "-", 
                           label = f'MAD = {round(signal_stats["MAD"], 3)}')
                ax.legend()

        for dist_idx, dist_name in enumerate(plot_distributions):
            # Get rows for the fitted distribution
            dist_rows = dist_stats_df[dist_stats_df["distribution"] == dist_name]
            # Find all parameter types (one per row)
            dist_param_types = np.intersect1d(np.unique(dist_rows["param_type"]), param_types)

            if dist_name in alt_dist_titles:
                dist_title = alt_dist_titles[dist_name]
            else:
                # Convert name to title case
                dist_title = dist_name.title()

            # Get the colour
            dist_colour = dist_colours[dist_name]

            for param_type in dist_param_types:
                # Select the axis to add the distribution to
                ax = axes[param_types_ax_idxs[param_type]]

                # Find stats for the distribution with specific parameters
                param_row = dist_rows[dist_rows["param_type"] == param_type]
                location = param_row["location"].iloc[0]
                scale = param_row["scale"].iloc[0]

                # Get the distribution's probabilty
                dist = self.scipy_distributions[dist_name]
                pdf = dist.pdf(x, loc = location, scale = scale)
                # Plot the fitted distribution
                ax.plot(x, pdf, label = dist_title, linewidth = 1.5, color = dist_colour, alpha = 0.9)

                if plot_zone_threshold:
                    # Translate from probability to threshold to compare against signal
                    zone_threshold = self.calculateZoneThreshold(sample = bw_idx,
                                                                 dist_name = dist_name,
                                                                 location = location,
                                                                 scale = scale,
                                                                 zone_probability = zone_probability,
                                                                 reverse_transform = False)
                    # Add the threshold line
                    ax.axvline(x = np.abs(zone_threshold), color = dist_colour,
                               alpha = line_alpha, linestyle = "--", 
                               label = f"{dist_title} Zone Threshold")

        # Set title above plots
        if title:
            plt.suptitle(title, fontweight = "bold")
        else:
            plt.suptitle(f"Fitted Distributions for {sample_name} ({chromosome})", 
                         fontweight = "bold")

        plt.legend()
        plt.tight_layout()

        if pdf_name:
            plt.savefig(os.path.join(self.output_directories["plots"], pdf_name), 
                        format = "pdf", bbox_inches = "tight")

        plt.show()

    def plotQQPlot(self, plot_sample, chromosome, param_type, signal = np.array([]), 
                   transform_signal = True, downsample = True, plot_distributions = [],
                   show_grid = True, title = "", plot_width = 15, plot_height = 3, pdf_name = ""):
        """
        Create quantile-quantile plot(s) per fitted distribution(s) for a specific parameter type.

        params:
            plot_sample:         Name of sample to plot values for.
            chromosome:          The chromosome to get plot values for.
            signal:              Optionally provide a signal to optionally transform and then create 
                                 the plot across.
            transform_signal:    Whether to transform the signal, e.g. by log transformation and/or 
                                 downsampling.
            downsample:          Set as True to downsample the signal.
            plot_distributions:  List of distributions to plot.
            param_type:          Distribution fitting parameter type to plot, e.g. "mean_fit" or "infer" 
                                 to use the type that had the overall lowest KS statistic across samples.
            show_grid:           Set as True to show a grid in the background.
            title:               Title to display at top of plot.
            plot_width:          Length of the plot.
            plot_height:         Height of the plot.
            pdf_name:            To save plot to PDF, set this as a file name.
        """

        custom_sample_names = self.getSampleNames(return_custom = True)

        if isinstance(plot_sample, str):
            bw_idx = self.sampleToIndex(plot_sample)
                
        elif isinstance(plot_sample, (int, np.integer)):
            # Use as sample ID
            bw_idx = plot_sample

            if bw_idx < 0 or bw_idx >= len(custom_sample_names):
                raise IndexError(f"Plot sample was set as a sample ID ({bw_idx}) that is out of range")
            
            plot_sample = custom_sample_names[bw_idx]

        if pdf_name:
            if not pdf_name.endswith(".pdf"):
                pdf_name = pdf_name + ".pdf"

            # Create directory to store plots
            os.makedirs(self.output_directories["plots"], exist_ok = True)

        plot_values = self.createDistributionPlotValues(bw_idx = bw_idx,
                                                        chromosome = chromosome,
                                                        signal = signal,
                                                        transform_signal = transform_signal,
                                                        downsample = downsample,
                                                        plot_distributions = plot_distributions,
                                                        calculate_stats = False)

        # Sanitised input distributions
        plot_distributions = plot_values["plot_distributions"]
        # The signal to plot a histogram of
        signal_to_fit = plot_values["signal_to_fit"]
        # Distribution statistics for the sample's chromosome
        dist_stats_df = plot_values["dist_stats_df"]
        del plot_values

        if param_type == "infer":
            best_zone_distribution = self.best_zone_distribution

            if best_zone_distribution is None:
                try:
                    self.inferBestDistribution()
                    best_zone_distribution = self.best_zone_distribution
                except Exception as e:
                    print(f"Tried to infer best distribution parameters, but failed due to exception: "
                          f"{e}")
                    return None

            param_type = list(best_zone_distribution.values())[0][0]

            if self.verbose > 0:
                print(f'Using inferred param_type "{param_type}"')

        else:
            # Check if parameter type is supported (ignoring suffix)
            param_type = param_type.lower().removesuffix("_fit")
            supported_param_types = [p.removesuffix("_fit") for p in self.param_types]
            if param_type in supported_param_types:
                all_param_types = [p.removesuffix("_fit") for p in np.unique(dist_stats_df["param_type"])]
                if param_type not in all_param_types:
                    raise ValueError(f'param_type "{param_type}" was not found in the distribution '
                                     f"statistics for {chromosome} {plot_sample}")
            else:
                raise ValueError(f'param_type "{param_type}" is not supported.\n'
                                f"Options include: {', '.join(self.param_types)}")

            # Add suffix back
            param_type += "_fit"

        # Create a subplot per distribution
        n_dists = len(plot_distributions)
        fig, axes = plt.subplots(1, n_dists, figsize = (plot_width, plot_height))
        
        if n_dists == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        alt_dist_titles = {"norm": "Normal",
                           "gumbel_l": "Gumbel Left",
                           "gumbel_r": "Gumbel Right"}

        for ax_idx, dist_name in enumerate(plot_distributions):
            # Get the axis to plot on
            ax = axes[ax_idx]

            # Get row for the fitted distribution with the specific parameters
            dist_row = dist_stats_df[(dist_stats_df["distribution"] == dist_name) &
                                     (dist_stats_df["param_type"] == param_type)]
            # Set the parameters
            location = dist_row["location"].iloc[0]
            scale = dist_row["scale"].iloc[0]
            params = (location, scale)

            # Create the Q-Q plot
            stats.probplot(signal_to_fit, dist = dist_name, sparams = params, plot = ax)

            if dist_name in alt_dist_titles:
                dist_title = alt_dist_titles[dist_name]
            else:
                # Convert name to title case
                dist_title = dist_name.title()

            ax.set_title(dist_title)
            ax.grid(show_grid)

            if ax_idx > 0:
                # Hide the y-axis label
                ax.tick_params(left = False, labelleft = False)
                ax.set_ylabel("")

        if title == "":
            title = f'Q-Q Plot{"s" if n_dists > 0 else ""} for {plot_sample} ({chromosome}) '
            title += f'- {param_type.title().replace("_", " ")}'
        
        plt.suptitle(title, fontweight = "bold")
        plt.tight_layout()

        if pdf_name:
            plt.savefig(os.path.join(self.output_directories["plots"], pdf_name), 
                        format = "pdf", bbox_inches = "tight")

        plt.show()