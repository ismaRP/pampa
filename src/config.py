#!/usr/bin/env python3

import json
import sys
import os

def parse_config_file():
    config_file = os.path.dirname(sys.argv[0]) + "/config.json"
    with open(config_file) as json_file:
        data = json.load(json_file)
        return data
# TODO: add PTM definition here
# TODO: add warning when the file is not found

def sort_headers(set_of_headers):
    config_file = os.path.dirname(sys.argv[0]) + "/config.json"
    with open(config_file) as json_file:
        data = json.load(json_file)
    list_of_selected_headers=[]
    list_of_other_headers=[]
    list_of_config_headers=[]
    for key in data:
        if  isinstance(data[key], list):
            list_of_config_headers.extend(data[key])
    for element in set_of_headers:
        if element in list_of_config_headers:
            list_of_selected_headers.append((list_of_config_headers.index(element), element))
        else:
            list_of_other_headers.append(element)
    list_of_selected_headers.sort(key=lambda x:x[0])
    return [z[1] for z in list_of_selected_headers]+list_of_other_headers
