#!/usr/bin/env python3

#########
# title: runphysicell.py
#
# language: python3
# license: BSD 3-Clause
# date: 2023-07-07
# author: bue
#
# description:
# python3 script to run physicell parameter scans.
# + https://github.com/MathCancer/PhysiCell
########


# libraries
import argparse
import json
import os
import re
import sys
import lxml.etree as ET
import pandas as pd


# functions
def runphysicell(s_settingxml, s_rulecsv, s_seedingcsv, s_paramjson, s_parammanipu, s_take):
    """
    input:
        check out argparse in __main__

    output:
        runs physicell.
        if manipulated, new PhysiCell_settings_run.xml
        if manipulated, new rules_run.csv

    description:
        code to run physicell, possibly with manipulated parameter settings.
    """
    print('runphysicell function ...')

    # handel input
    s_runsettingxml = s_settingxml
    if (s_parammanipu is None):
        s_pmanipu = 'none'
    else:
        s_pmanipu = re.sub('[^A-Za-z0-9.]', '', s_parammanipu.replace('-','N').replace('+','P'))

    # initialize flags
    b_setting = False
    b_rule = False

    # load setting xml file
    x_tree = ET.parse(s_settingxml)
    x_root = x_tree.getroot()

    # manipulate seeding csv xml element
    s_xpath_seed = './/initial_conditions/cell_positions//filename'
    lx_element = x_root.findall(s_xpath_seed)
    if len(lx_element) != 1:
        sys.exit(f'Error @ runphysicell : in {s_settingxml}, found no or more than one element {lx_element} for xpath {s_xpath_seed}.')
    if (s_seedingcsv is None):
        # get sedding csv filename
        s_seedingcsv = lx_element[0].text
    else:
        # manipulate sedding csv filename element
        b_setting = True
        lx_element[0].text = s_seedingcsv
    s_seeding = s_seedingcsv.replace('\\','/').split('/')[-1].replace('.csv','')

    # manipulate output folder xml element
    s_xpath_output = './/save/folder'
    lx_element = x_root.findall(s_xpath_output)
    if len(lx_element) != 1:
        sys.exit(f'Error @ runphysicell : in {s_settingxml}, found no or more than one element {lx_element} for xpath {s_xpath_output}.')
    s_out = f'output_{s_seeding}_{s_pmanipu}_{s_take}'
    if (lx_element[0].text != s_out):
        b_setting = True
        lx_element[0].text = s_out

    # manipulare xml and csv parameters
    if b_setting or not (s_rulecsv is None) or not (s_parammanipu is None):

        # load rule csv file
        s_xpath_rule = './/cell_rules/rulesets//filename'
        lx_element = x_root.findall(s_xpath_rule)
        if len(lx_element) != 1:
            sys.exit(f'Error @ runphysicell : in {s_settingxml}, found no or more than one element {lx_element} for xpath {s_xpath_rule}.')
        if (s_rulecsv is None):
            # get rule csv filename
            s_rulecsv = lx_element[0].text
        else:
            # manipulate rule csv filename element
            b_setting = True
            lx_element[0].text = s_rulecsv
        df_rule = pd.read_csv(s_rulecsv, comment='/', header=None)
        df_rule.columns = ['celltype','signal','direction','behavoir','saturation','halfmax','hillpower','applytodead']
        df_rule.index = [re.sub('[^A-Za-z0-9]','',s) for s in (df_rule.celltype + df_rule.signal + df_rule.direction + df_rule.behavoir)]

        # manipulate parameters in setting xml and rule csv
        if not (s_parammanipu is None):
            s_pmanipu = re.sub('[^A-Za-z0-9.]', '', s_parammanipu.replace('-','N').replace('+','P'))

            # load parameter json file
            d_param = json.load(open(s_paramjson))
            # the loop
            for l_param in d_param[s_parammanipu]:
                # manipulate parameter setting xml element
                if (l_param[0] == 'xml'):
                    lx_element = x_root.xpath(l_param[1])  # xpath
                    for x_element in lx_element:
                        x_element.text = str(l_param[2])  # value
                # manipulate parameter rule csv element
                elif (l_param[0] == 'csv'):
                    df_rule.loc[l_param[1],l_param[2]] = l_param[3]  # row, column, value
                # error handling
                else:
                    sys.exit(f'Error @ runphysicell : in s_paramjson file {s_paramjson} for s_parammanipu label {s_parammanipu} unknown manipulation file type detected "{l_param[0]}". known are csv and xml.')

        # write rule csv file
        if b_rule:
            s_runrulecsv = f"{s_settingxml.replace('.xml','_')}{s_seeding}_{s_pmanipu}.xml"
            df_rule.to_csv(s_runrulecsv, index=False, header=False)

        # write settings xml file
        if b_setting:
            s_runsettingxml = f"{s_settingxml.replace('.xml','_')}{s_seeding}_{s_pmanipu}.xml"
            x_tree.write(
                s_runsettingxml,
                xml_declaration='<?xml version="1.0" encoding="UTF-8"?>'
            )

    # run dmc
    os.mkdir(s_out)
    os.system(f'physicell {s_runsettingxml}')
    print('runphysicell finished!')


# main
if __name__ == '__main__':
    print('runphysicell script ...')
    # argv
    parser = argparse.ArgumentParser(
        prog = 'runphysicell',
        description = 'manipulates template seetings.xml and rules.csv files for physicell parameter scan.',
        epilog = 'may the force be with you!',
    )
    # model
    parser.add_argument(
        'settingxml',
        #nargs = 1,
        help = 'PhysiCell_seetings.xml model file.'
    )
    parser.add_argument(
        '-r', '--rule',
        default = None,  # default specified in settings.xml
        help = 'to specify an alternative rules.csv model file to the one specified in the settings.xml file.',
        # rules.csv header 1.13.1: celltype,signal,direction,behavoir,saturation,halfmax,hillpower,applytodead
    )
    # parameter manipulation
    parser.add_argument(
        'manipu',
        nargs = '?',
        default = None,  # no parameter manipulation
        help = 'to specify which parameter settings set in parameter.json to use to manipulate the template model seetings.xml and rules.csv files.'
    )
    # seeding
    parser.add_argument(
        '-i', '--init',
        default = None,  # default specified in settings.xml
        help = 'to specify an alternative cells.csv seeding (initial condition) file to the one specified in the settings.xml file.',
    )
    # parameters
    parser.add_argument(
        '-p', '--param',
        default = None,
        help = 'to specify a parameter scan json file.'
    )
    # takes
    parser.add_argument(
        '-t', '--take',
        #type = int,
        default = 'a',
        help = 'to specify the take per parameter setting label.',
    )
    # parse arguments
    args = parser.parse_args()
    print(args)

    # parameter sanity check
    if not (args.manipu is None):
        if (args.param is None):
            sys.exit(f'Error @ runphysicell : parameter set {args.manipu} for manipulation specified, but -p path/to/parameter.json argumnet, where the set can be found, is missing!')

    # processing
    runphysicell(
        s_settingxml = args.settingxml,
        s_rulecsv = args.rule,
        s_seedingcsv = args.init,
        s_paramjson = args.param,
        s_parammanipu = args.manipu,
        s_take = args.take,
    )

