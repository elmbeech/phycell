#!/usr/bin/env python3

# library
import argparse
import json
import os
import sys
import lxml.etree as ET
import pandas as pd


# functions

#s_parameter = '$s_parameter'
#s_xpath_parameter = '$xp_parameter'
#r_parameter = '$r_parameter'

def runphysicell(s_settingxml, s_rulecsv, s_seedingcsv, s_paramjson, s_parammanipu, i_take):
    """
    input: 
        check out argparse in __main__

    output:

    description:
    """
    print("yepa!")

    # processing
    s_runsettingxml = s_settingxml
    if not(s_rulecsv is None) or not(s_seedingcsv is None) or not(s_parammanipu is None):

        # load parameter json file


        # load setting xml file
        x_tree = ET.parse(s_settingxml)
        x_root = x_tree.getroot()

        # maniplate rules xml element
        s_xpath_rule = './/cell_rules/rulesets//filename'
        lx_element = x_root.findall(s_xpath_rule)
        if len(lx_element) != 1:
            print('xpath:', s_xpath_rule)
            sys.exit('Error @ runPhysiCell : in {s_settingxml}, found non or more than one element for xpath.')
        lx_element[0].text = s_rulecsv

        # load rule csv file


        # manipulate initial condition xml element
        s_xpath_init = './/initial_conditions/cell_positions//filename'
        lx_element = x_root.findall(s_xpath_init)
        if len(lx_element) != 1:
             print('xpath:', s_xpath_init)
             sys.exit('Error @ runPhysiCell : in {s_settingxml}, found non or more than one element for xpath.')
        lx_element[0].text = s_seedingcsv

        # manipulate output folder xml element
        s_xpath_output = './/save/folder'
        lx_element = x_root.findall(s_xpath_output)
        if len(lx_element) != 1:
            print('xpath:', s_xpath_output)
            sys.exit(f'Error @ runPhysiCell : in {s_settingxml}, found non or more than one element for xpath.')
        s_out = f"ouput_{s_seedingcsv.replace('.csv','')}_{s_parameter}{r_parameter}_{i_take}"
        os.mkdir(s_out)
        lx_element[0].text = s_out

        # manipulate parameter setting xml element
        if s_xpath_parameter != 'None':
            lx_element = x_root.xpath(s_xpath_parameter)
            for x_element in lx_element:
                x_element.text = r_parameter

        # manipulate parameter rule csv element
        # process rule.csv file


        # write rule csv file
        #s_nfsetting = f"{s_settingxml.replace('.xml','_')}{s_seedingcsv.replace('.csv','')}_{s_parameter}{r_parameter}.xml"
        #x_tree.write(
        #    s_nfsetting,
        #    xml_declaration='<?xml version="1.0" encoding="UTF-8"?>'
        #)

        # write settings xml file
        #s_runsettingxml = f"{s_settingxml.replace('.xml','_')}{s_seedingcsv.replace('.csv','')}_{s_parameter}{r_parameter}.xml"
        #x_tree.write(
        #    s_runsettingxml,
        #    xml_declaration='<?xml version="1.0" encoding="UTF-8"?>'
        #)

    # run dmc
    #os.system(f'physicell {s_runsettingxml}')


# main
if __name__ == '__main__':
    # argv
    parser = argparse.ArgumentParser(
        prog = 'runphysicell',
        description = 'manipulates templare seetings.xml and rules.csv file for physicell parameter scan',
        epilog = 'may the force be with you!',
    )
    # model
    parser.add_argument(
        'settingxml', 
        nargs = 1,
        help = 'PhysiCell_seetings.xml model file.'
    )
    parser.add_argument(
        '-r', '--rule', 
        default = None,  # default specified in settings.xml
        help = 'specify rules.csv model file',
        # rules.csv header 1.13.1: celltype,signal,direction,behavoir,saturation,halfmax,hillpower,applytodead
    )
    # parameter manipulation
    parser.add_argument(
        'manipu',
        nargs = '?',
        default = None,  # no parameter manipulation
        help = 'specify which parameter settings set from the parameter.json to user to manipulate the template model seetings.xml and rules.csv file.'
    )
    # seeding
    parser.add_argument(
        '-i', '--init', 
        default = None,  # default specified in settings.xml
        help = 'specify cells.csv seeding file',
    )
    # parameters
    parser.add_argument(
        '-p', '--param',
        default = 'parameter.json',
        help = 'specify parameter scan json file.'
    )
    # takes
    parser.add_argument(
        '-t', '--take',
        type = int,
        default = 1,
        help = 'specify how many takes to run per parameter setting.',
    )
    # parse arguments
    args = parser.parse_args()
    print(args)

    # parameter sanity check
    if not (args.manipu is None):
        print("args.manipu:", args.manipu)
        if not os.path.exists(args.param):
            sys.exit(f"Error @ runphysicell : {args.param} file not found. utilize the -p or --param argument to expliciely specify path/to/parameter.json!")

    # processing
    runphysicell(
        s_settingxml = args.settingxml, 
        s_rulecsv = args.rule, 
        s_seedingcsv = args.init,
        s_paramjson = args.param,
        s_parammanipu = args.manipu,
        i_take = args.take,
    )

