#create object document for extracting the ms1 related information
class Document:

    def __init__(self, peak_id, mz, rt, intensity, charge, ms2_peak_list):

        self.peak_id = peak_id
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.charge = charge

        self.spectra = ms2_peak_list


#create object peak for extracting ms2 related information

class Peak:

    def __init__(self, peak_id, parent_id, mz, rt, intensity, charge):
        self.peak_id=peak_id
        self.parent_id = parent_id
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.charge = charge
