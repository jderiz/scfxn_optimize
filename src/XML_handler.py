import copy
import re


class XML_Handler():

    def __init__(self, xml_script: str, tags: list):
        self.input = xml_script
        self.tag_list = tags

    def get_script_for_config(config_list: list) -> str:

        for tag, new_value in zip(self.tag_list, config_list):
            # substitute the new value from optimizer in xml script
            modified = re.sub(r''+tag+'='+'\".*\"', tag +
                              '='+str(new_value), self.input)

        return modified
