import copy
import re

import pyrosetta.rosetta.protocols.rosetta_scripts as prs_scripts
import pyrosetta.rosetta.utility.options as options


class XML_Handler():

    def __init__(self, xml_script: str, tag_list: list, option_file_path: str):
        self.input = xml_script
        self.xml_for_config = "/tmp/pyrosetta/xml_for_config.xml"
        self.tag_list = tag_list
        self.option_file_path = option_file_path

    def get_script_for_config(config_list: list) -> str:

        for tag, new_value in zip(self.tag_list, config_list):
            # substitute the new value from optimizer in xml script
            modified = re.sub(r''+tag+'='+'\".*\"', tag +
                              '='+str(new_value), self.input)
            # write to disk
            with open(self.xml_for_config, 'w') as handle:
                handle.write(modified)

        return modified

    def objective(config, args):
        option_collection = prs.rosetta.basic.options.initialize()  # TODO: figure out
        option_collection.load_options_from_file(self.option_file_path)

        parser = prs.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()
        parsed_protocol: prs.rosetta.protocols.rosetta_scripts.ParsedProtocol = parser.generate_mover_xml_string(
            self.get_script_for_config(config), option_collection)
        # load modified xml file from disk
        parsed_protocol: prs.rosetta.protocols.rosetta_scripts.ParsedProtocol = parser.generate_mover(
            xml_fname=self.xml_for_config,
            options=self.option_file_path,
            script_vars=self.script_vars_vector
        )
        parsed_protocol.apply(pose)
