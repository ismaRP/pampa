
import os
import pyopenms as oms
import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import product
import pandas as pd
import textwrap
import seaborn as sns
import itertools as it

# local import
from src import message

def set_filter_parameters(filter, params):
    param = oms.Param()
    for key, value in params.items():
        param.setValue(key, value)
    filter.setParameters(param)


def fix_peaks_spectrum(peaks_spec, profile_spec):
    mz_peaks, _ = peaks_spec.get_peaks()
    _, int_prof = profile_spec.get_peaks()
    peak_indices = [profile_spec.findHighestInWindow(m, 0.1, 0.1) for m in mz_peaks]
    int_peaks = int_prof[peak_indices]
    peaks_spec.set_peaks((mz_peaks, int_peaks))
    return peaks_spec


def fix_peaks_int(peaks_exp, profile_exp):
    assert peaks_exp.getNrSpectra() == profile_exp.getNrSpectra()
    n = peaks_exp.getNrSpectra()
    new_spectra = []
    for i in range(n):
        prof_spec = profile_exp.getSpectrum(i)
        peaks_spec = peaks_exp.getSpectrum(i)
        peaks_spec = fix_peaks_spectrum(peaks_spec, prof_spec)
        new_spectra.append(peaks_spec)
    peaks_exp.setSpectra(new_spectra)
    return peaks_exp

def map_repl_spectra_id(deiso_exp):
    key_spectra = dict()
    for spec in deiso_exp.getSpectra():
        spec_name = spec.getNativeID()
        sample_name = spec_name.split('_')[0]
        key_spectra[spec_name] = sample_name
    return key_spectra


def map_repl_file(repl_to_sample):
    repl_to_sample = csv.reader(open(repl_to_sample, "r"))
    key_spectra = dict()
    for row in repl_to_sample:
        key_spectra[row[0]] = row[1]
    return key_spectra

class PreprocessingConsumer:
    def __init__(self, consumer_deiso, params, consumer_iso=None, normalize=False):
        self._internal_consumer = consumer_deiso
        self._internal_consumer_iso = consumer_iso
        self.sg_filter = oms.SavitzkyGolayFilter()
        self.morph_filter = oms.MorphologicalFilter()
        self.normalizer = oms.Normalizer()
        self.peak_picker = oms.PeakPickerIterative()
        self.normalize = normalize
        # Set parameters for each filter based on params dictionary
        for f, par in params.items():
            if f == 'SGF':
                set_filter_parameters(self.sg_filter, par)
            if f == 'BC':
                set_filter_parameters(self.morph_filter, par)
            if f == 'PP':
                set_filter_parameters(self.peak_picker, par)
            if f == 'NORM':
                set_filter_parameters(self.normalizer, par)

    def setExperimentalSettings(self, s):
        self._internal_consumer.setExperimentalSettings(s)

    def setExpectedSize(self, a, b):
        self._internal_consumer.setExpectedSize(a, b)

    def consumeChromatogram(self, c):
        self._internal_consumer.consumeChromatogram(c)

    def consumeSpectrum(self, s):
        # Apply the Savitzky-Golay filter
        self.sg_filter.filter(s)
        # Apply the morphological filter
        self.morph_filter.filter(s)
        # Apply peak picking
        picked_spectrum = oms.MSSpectrum()
        self.peak_picker.pick(s, picked_spectrum)
        picked_spectrum = fix_peaks_spectrum(picked_spectrum, s)
        if self._internal_consumer_iso is not None:
            self._internal_consumer_iso.consumeSpectrum(picked_spectrum)
        oms.Deisotoper.deisotopeAndSingleCharge(
            spectra=picked_spectrum,
            fragment_tolerance=50 , fragment_unit_ppm=True,
            min_charge=1, max_charge=1, keep_only_deisotoped=True,
            min_isopeaks=3 , max_isopeaks=6, make_single_charged=False,
            annotate_charge=False, annotate_iso_peak_count=False,
            use_decreasing_model=True, start_intensity_check=2 , add_up_intensity=False
        )
        if self.normalize:
            self.normalizer.filterPeakSpectrum(picked_spectrum)
        # Pass the processed spectrum to the internal consumer
        self._internal_consumer.consumeSpectrum(picked_spectrum)


def merge_replicates(deiso_exp, map_repl):
    # Start the merger
    repl_merger = oms.SpectraMerger()
    # Set the parameters
    merger_params = {
        # Tolerance for binning peaks between replicates
        'mz_binning_width': 0.4,
        'mz_binning_width_unit': 'Da',
        # Here we are telling to merge blocks of up to 5 spectra, which should be enough
        # since we are working with blocks of 3 replicates
        'block_method:rt_block_size': 5,
        # This parameter set a maximum RT range for a block
        # Using 0.0 will disable any RT restriction, since we're not even working with RTs
        'block_method:rt_max_length': 0.0,
        # Here we tell the algorithm to only go through MS1 spectra for merging,
        # which is what we have.
        'block_method:ms_levels': [1]
    }
    # Pass parameters to the merger
    set_filter_parameters(repl_merger, merger_params)

    merged_spectra = []
    for k, spec_gr in it.groupby(deiso_exp.getSpectra(), lambda x: map_repl.get(x.getNativeID())):
        spec_gr = list(spec_gr)
        # If all 3 replicates are empty, create an empty merged spectrum manually
        if all([s.size() == 0 for s in spec_gr]):
            merged_spectrum = oms.MSSpectrum()
        else:
            # Create a temp MSExperiment with the replicates
            sample_exp = oms.MSExperiment()
            sample_exp.setSpectra(spec_gr)
            # Perform merging
            repl_merger.mergeSpectraBlockWise(sample_exp)
            merged_spectrum = sample_exp.getSpectrum(0)
        merged_spectrum.setNativeID(k)
        merged_spectrum.setName(k)
        merged_spectra.append(merged_spectrum)
    # Initialize a new merged MSExperiment
    allsamples_exp = oms.MSExperiment()
    # Add the list of merged spectra to it
    allsamples_exp.setSpectra(merged_spectra)
    return allsamples_exp

def to_peaklist(spectra, prep_params, repl_to_sample, output, make_plot=False):
    output_dir, output_file = os.path.split(output)
    deiso_writer = oms.PlainMSDataWritingConsumer(os.path.join(output_dir, '/proc_deisotoped.mzML'))
    nodeiso_writer = None
    if make_plot:
        nodeiso_writer = oms.PlainMSDataWritingConsumer(output_dir + '/proc_nodeisotoped.mzML')
    prep_consumer = PreprocessingConsumer(
        deiso_writer, prep_params['params'], nodeiso_writer, normalize=False)

    if os.path.isfile(spectra):
        spectra_list = [spectra]
    elif os.path.isdir(spectra):
        spectra_list = [os.path.join(spectra, f) for f in os.listdir(spectra) if f.endswith('mzML')]
    for m in spectra_list:
        oms.MzMLFile().transform(m, prep_consumer)

    del deiso_writer, nodeiso_writer, prep_consumer

    deiso_exp = oms.MSExperiment()
    oms.MzMLFile().load(os.path.join(output_dir, '/proc_deisotoped.mzML'), deiso_exp)
    deiso_exp.clearMetaDataArrays()
    if make_plot:
        nodeiso_exp = oms.MSExperiment()
        oms.MzMLFile().load(output_dir + '/proc_nodeisotoped.mzML', nodeiso_exp)
        nodeiso_exp.clearMetaDataArrays()
        plot_npeaks(deiso_exp, nodeiso_exp, output_dir)

    # Merge replicates
    if repl_to_sample is not None:
        map_repl = map_repl_file(repl_to_sample)
    else:
        map_repl = map_repl_spectra_id(deiso_exp)
    merged_spectra = merge_replicates(deiso_exp, map_repl)
    merged_writer = oms.PlainMSDataWritingConsumer(output)
    merged_writer.setExperimentalSettings(merged_spectra.getExperimentalSettings())
    for s in merged_spectra.getSpectra():
        merged_writer.consumeSpectrum(s)
    del merged_writer


def abline(slope, intercept, color, ax=None, label=None):
    if ax is None:
        ax = plt.gca()
    x_vals = np.array(ax.get_xlim())
    y_vals = intercept + slope * x_vals
    ax.plot(x_vals, y_vals, '--', color=color, label=label)


def plot_npeaks(deiso_exp, nodeiso_exp, output_dir):
    npeaks_deiso = []
    npeaks_nodeiso = []
    for deiso_s, nodeiso_s in zip(deiso_exp.getSpectra(), nodeiso_exp.getSpectra()):
        npeaks_deiso.append(deiso_s.size())
        npeaks_nodeiso.append(nodeiso_s.size())
    fig, ax = plt.subplots(
        nrows=1, ncols=1, figsize=(6, 4), layout='constrained')
    ax.scatter(npeaks_nodeiso, npeaks_deiso, alpha=0.7)
    abline(1 / 5, 0, '#1f77b4', ax=ax, label='5 to 1')
    abline(1 / 6, 0, '#ff7f0e', ax=ax, label='6 to 1')
    ax.set_xlabel('# of peaks')
    ax.set_ylabel('# of monoisotopic peaks')
    ax.legend(loc='best')
    fig.savefig(output_dir + '/npeaks.png')


def main(spectra, output, prep_params, repl_to_sample, make_plot):

    output_dir, output_file = os.path.split(output)
    if len(output_dir)>0:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = os.getcwd()

    message.configure(output_dir)
    to_peaklist(spectra, prep_params, repl_to_sample, output, make_plot)

if __name__ == '__main__':
    main()
    

