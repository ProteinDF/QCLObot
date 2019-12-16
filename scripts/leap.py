#!/usr/bin/env python
# -*- coding: utf-8 -*-

class QcLeap(object):
    def __init__(self):
        self._leaprcs = [
            'leaprc.ff03.r1',
            'leaprc.gaff']
        self._pdb_models = {}

    # ------------------------------------------------------------------
    def set_models(self, model_name, pdb_file):
        self._pdb_models.set(model_name, pdb_file)

    
        
    # ------------------------------------------------------------------
    def _output_leapin(self):
        leapin_templ = """
        proteinBox = copy protein
        
        quit
        """

        leapin_contents = ""
        leapin_contents += self._prepare_leapin_source()
        leapin_contents += self._prepare_leapin_pdbs()
        leapin_contents += self._prepare_leapin_ssbonds()
        leapin_contents += self._prepare_leapin_save_params()
        leapin_contents += self._prepare_leapin_save_models()
        
    def _prepare_leapin_source(self):
        output = ""
        for leaprc in self._leaprcs:
            output += "source {}\n".format(leaprc)
        return output
        
    def _prepare_leapin_pdbs(self):
        output = ""
        for model_name, pdb_file in self._pdb_models:
            output += "{} = loadPdb {}\n".format(model_name, pdb_file)

        return output

    def _prepare_leapin_ssbonds(self):
        output = ""
        for ssbond in self._ssbonds:
            output += "\n"
        return output

    def _prepare_leapin_solvatecap(self, model):
        output = ""
        for solvatecap in self._solvatecaps:
            output += "solvateCap {model} TIP3PBOX ${solcap_center} ${solcap_closeness}\n".format(
                model = solvatecap["model"],
                wat_model = solvatecap["wat_model"], # TIP3PBOX etc.
                solcap_center = solvatecap["solcap_center"],
                solcap_closeness = solvatecap["solcap_closeness"]
            )

        return output
        
    def _prepare_leapin_save_params(self):
        output = ""
        output += "saveAmberParm proteinBox md1.prmtop md1.inpcrd\n".format(self.prmtop_file, self.inpcrd)
        return output

    def _prepare_leapin_save_models(self):
        output = ""
        for model_name, pdb_file in self._save_models:
            output += "savepdb {model_name} {pdb_file}\n".format(
                model_name = model_name,
                pdb_file = pdb_file)
        return output


def main():
    leap = QcLeap()
    leap.load_pdb("sample.pdb", sample)
