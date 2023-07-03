#########
# title: pyMCDS.py
#
# language: python3
# date: 2022-08-22
# license: BSD-3-Clause
# authors: Patrick Wall, Randy Heiland, Paul Macklin, Elmar Bucher
#
# description:
#     pyMCDS.py definds an object class, able to load and access
#     within python a single time step from the PhysiCell model output folder.
#     pyMCDS.py was forked from the original PhysiCell-Tools python-loader
#     implementation and further developed.
#########


# load library
import pandas as pd
import xml.etree.ElementTree as ET


# const physicell codec
# implemation based on PhysiCell/core/PhysiCell_constants.h.
ds_cycle_model = {
    '0' : 'advanced_Ki67_cycle_model',
    '1' : 'basic_Ki67_cycle_model',
    '2' : 'flow_cytometry_cycle_model',
    '3' : 'live_apoptotic_cycle_model',
    '4' : 'total_cells_cycle_model',
    '5' : 'live_cells_cycle_model',
    '6' : 'flow_cytometry_separated_cycle_model',
    '7' : 'cycling_quiescent_model',
}
ds_death_model = {
    '100' : 'apoptosis_death_model',
    '101' : 'necrosis_death_model',
    '102' : 'autophagy_death_model',
    '9999' : 'custom_cycle_model',
}

ds_cycle_phase = {
    '0' : 'Ki67_positive_premitotic',
    '1' : 'Ki67_positive_postmitotic',
    '2' : 'Ki67_positive',
    '3' : 'Ki67_negative',
    '4' : 'G0G1_phase',
    '5' : 'G0_phase',
    '6' : 'G1_phase',
    '7' : 'G1a_phase',
    '8' : 'G1b_phase',
    '9' : 'G1c_phase',
    '10' : 'S_phase',
    '11' : 'G2M_phase',
    '12' : 'G2_phase',
    '13' : 'M_phase',
    '14' : 'live',
    '15' : 'G1pm_phase',
    '16' : 'G1ps_phase',
    '17' : 'cycling',
    '18' : 'quiescent',
    '9999' : 'custom_phase',
}
ds_death_phase = {
    '100' : 'apoptotic',
    '101' : 'necrotic_swelling',
    '102' : 'necrotic_lysed',
    '103' : 'necrotic',
    '104' : 'debris',
}


# object classes
class Parameter:
    """
    input:
        s_path: string; default '.'
            relative or absolute path to the directory where
            the PhysiCell output files are stored.
    output:
        mcds: pyMCDS class instance
            all fetched content is stored at mcds.data.

    description:
        pyMCDS.__init__ will generate a class instance with a
        dictionary of dictionaries data structure that contains all
        output from a single PhysiCell model time step. furthermore,
        this class, and as such it's instances, offers functions
        to access the stored data.
        the code assumes that all related output files are stored in
        the same directory. data is loaded by reading the xml file
        for a particular time step and the therein referenced files.
    """
    def __init__(self, s_path='.'):
        self.data = self._read_setting_xml(s_path)


    ## PARAMETER RELATED FUNCTIONS ##

    def get_parameter_dict(self):
        """
        input:
            self: pyMCDS class instance.

        output:
            d_parameter: dictionary
            dictionary, mapping each tracked paramater to its input value.

        description:
            function retunes a dictionary that maps
            the models input parameter and values.

        attention:
        this is an unofficial function and can with any version be dropped.
        this will for sure happen, when cell_type and substrate id label
        mapping dictionaries are provided in the regular output xml,
        to keep the code lean and agile!
        """
        d_parameter = self.data['setting']['parameters'].copy()
        return(d_parameter)

    def get_unit_se(self):
        """
        input:
            self: pyMCDS class instance.

        output:
            se_unit: pandas series
            series lists all tracked parameter variables that have units,
            and maps them to their unit.

        description:
            function returns a series that lists all tracked
            parameter variables that have units and their units.

        attention:
        this is an unofficial function and can with any version be dropped.
        this will for sure happen, when cell_type and substrate id label
        mapping dictionaries are provided in the regular output xml,
        to keep the code lean and agile!
        """
        se_unit = pd.Series(self.data['setting']['units'])
        se_unit.index.name = 'parameter'
        se_unit.name = 'unit'
        se_unit.sort_index(inplace=True)
        return(se_unit)

    def get_rule_df(self):
        """
        input:
            self: pyMCDS class instance.

        output:
            df_rules: pandas dataframe
            rules.csv loaded as datafarme.

        description:
            function returns the rule csv as a datafram.

        attention:
        this is an unofficial function and can with any version be dropped.
        this will for sure happen, when cell_type and substrate id label
        mapping dictionaries are provided in the regular output xml,
        to keep the code lean and agile!
        """
        df_rule = self.data['setting']['rules']
        return(df_rule)


    ## LOAD DATA ##

    def _read_setting_xml(self, s_path='.'):
        """
        input:
            self: pyMCDS class instance.

            s_path: string; default '.'
                relative or absolute path to the directory where
                the PhysiCell output files, are stored.
                especially the PhysiCell_setting.xml file.

        output:
            self: Parameter class instance with loaded data.

            d_mcds updated with:
            d_mcds['metadata']['cell_type'] dictionary
            d_mcds['metadata']['substrate'] dictionary
            d_mcds['setting']['parameters'] dictionary
            d_mcds['setting']['units'] dictionary
            d_mcds['setting']['rules'] pandas dataframe

        description:
            this internal function extract information form the
            PhysiCell_setting.xml file.
        """
        # generate data storage dictionary
        d_mcds = {}
        ds_substrate = {}
        ds_celltype = {}
        d_parameter = {}
        ds_unit = {}
        df_ruleset = None

        ### load Physicell_settings xml file ###
        s_xmlpathfile = s_path + '/PhysiCell_settings.xml'
        print(f'reading: {s_xmlpathfile}')
        x_tree = ET.parse(s_xmlpathfile)
        x_root = x_tree.getroot()

        ### find the overall node ###
        x_overall = x_root.find('overall')
        if not (x_overall is None):
            # dt_diffusion
            try:
                d_parameter.update({'dt_diffusion': float(x_overall.find('dt_diffusion').text)})
                ds_unit.update({'dt_diffusion': x_overall.find('dt_diffusion').get('units')})
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <overall><dt_diffusion> node missing.')
            # dt_mechanics
            try:
                d_parameter.update({'dt_mechanics': float(x_overall.find('dt_mechanics').text)})
                ds_unit.update({'dt_mechanics': x_overall.find('dt_mechanics').get('units')})
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <overall><dt_mechanics> node missing.')

            # dt_phenotype
            try:
                d_parameter.update({'dt_phenotype': float(x_overall.find('dt_phenotype').text)})
                ds_unit.update({'dt_phenotype': x_overall.find('dt_phenotype').get('units')})
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <overall><dt_phenotype> node missing.')
        else:
            print(f'Warning @ pyMCDS._read_setting_xml : <overall> node missing.')

        ### find the options node ###
        x_options = x_root.find('options')
        if not (x_options is None):
            #
            try:
                d_parameter.update({'legacy_random_points_on_sphere_in_divide': str(x_overall.find('legacy_random_points_on_sphere_in_divide').text).lower == 'true'})
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <options><legacy_random_points_on_sphere_in_divide> node missing.')
            #
            try:
                d_parameter.update({'virtual_wall_at_domain_edge': str(x_overall.find('virtual_wall_at_domain_edge').text).lower == 'true'})
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <options><virtual_wall_at_domain_edge> node missing.')
            #
            try:
                d_parameter.update({'disable_automated_spring_adhesions': str(x_overall.find('disable_automated_spring_adhesions').text).lower == 'true'})
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <options><disable_automated_spring_adhesions> node missing.')
        else:
            print(f'Warning @ pyMCDS._read_setting_xml : <options> node missing.')

        ### find the microenvironment node ###
        x_microenvironment = x_root.find('microenvironment_setup')
        # substrate loop
        for x_variable in x_microenvironment.findall('variable'):
            # basics
            i_id = int(x_variable.get('ID'))
            s_substrate = x_variable.get('name').replace(' ', '_')
            ds_substrate.update({i_id : s_substrate})

            # microenvironment physics
            x_microenvironment_physical = x_variable.find('physical_parameter_set')
            if not (x_microenvironment_physical is None):
                # diffusion_coefficient
                try:
                    d_parameter.update({f'{s_substrate}_diffusion_coefficient': float(x_microenvironment_physical.find('diffusion_coefficient').text)})
                    ds_unit.update({f'{s_substrate}_diffusion_coefficient': x_microenvironment_physical.find('diffusion_coefficient').get('units')})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{i_id}"><physical_parameter_set><diffusion_coefficient> node missing.')
                # decay_rate
                try:
                    d_parameter.update({f'{s_substrate}_decay_rate': float(x_microenvironment_physical.find('decay_rate').text)})
                    ds_unit.update({f'{s_substrate}_decay_rate': x_microenvironment_physical.find('decay_rate').get('units')})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{i_id}"><physical_parameter_set><decay_rate> node missing.')
            else:
                print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{i_id}"><physical_parameter_set> node missing.')
            # microenv initial
            try:
                d_parameter.update({f'{s_substrate}_initial_condition': float(x_variable.find('initial_condition').text)})
                ds_unit.update({f'{s_substrate}_initial_condition': x_variable.find('initial_condition').get('units')})
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{i_id}"><initial_condition> node missing.')
            # microenv dirichlet
            try:
                d_parameter.update({f'{s_substrate}_dirichlet_boundary_condition_enabled': str(x_variable.find('Dirichlet_boundary_condition').get('enabled')).lower == 'true'})
                d_parameter.update({f'{s_substrate}_dirichlet_boundary_condition': float(x_variable.find('Dirichlet_boundary_condition').text)})
                ds_unit.update({f'{s_substrate}_dirichlet_boundary_condition': x_variable.find('Dirichlet_boundary_condition').get('units')})
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{i_id}"><Dirichlet_boundary_condition> node missing.')
            x_microenvironment_dirichlet = x_variable.find('Dirichlet_options')
            if not (x_microenvironment_dirichlet is None):
                for x_dirichlet in x_microenvironment_dirichlet:
                    s_coor = x_dirichlet.get('ID')
                    d_parameter.update({f'{s_substrate}_dirichlet_boundary_value_{s_coor}_enabled': str(x_dirichlet.text).lower == 'true'})
                    d_parameter.update({f'{s_substrate}_dirichlet_boundary_value_{s_coor}': float(x_dirichlet.text)})
            else:
                print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{i_id}"><Dirichlet_options> node missing.')

        ### find the cell definition node ###
        # cell loop
        es_customdata = set()
        for x_celltype in x_root.find('cell_definitions').findall('cell_definition'):
            # basics
            i_id = int(x_celltype.get('ID'))
            s_celltype = x_celltype.get('name').replace(' ', '_')
            ds_celltype.update({i_id : s_celltype})

            # search for phenotype
            x_phenotype = x_celltype.find('phenotype')
            if not (x_phenotype is None):

                # phenotype live cycle model
                x_cycle = x_phenotype.find('cycle')
                if not (x_cycle is None):
                    d_parameter.update({f'{s_celltype}_cycle_model': ds_cycle_model[str(x_cycle.get('code'))]})
                    # duration
                    try:
                        for x_phase in x_cycle.find('phase_durations').findall('duration'):
                            s_index = ds_cycle_phase[str(x_phase.get('index'))]
                            d_parameter.update({f'{s_celltype}_cycle_phase_duration_{s_index}_fixed': str(x_phase.get('fixed_duration')).lower() == 'true'})
                            d_parameter.update({f'{s_celltype}_cycle_phase_duration_{s_index}': float(x_phase.text)})
                            ds_unit.update({f'{s_celltype}_cycle_phase_duration_{s_index}': x_phase.get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cycle><phase_durations> node missing.')
                    # rate
                    try:
                        for x_phase in x_cycle.find('phase_transition_rates').findall('rate'):
                            s_index = ds_cycle_phase[str(x_phase.get('start_index'))] + '_' + ds_cycle_phase[str(x_phase.get('end_index'))]
                            d_parameter.update({f'{s_celltype}_cycle_phase_transition_rate_{s_index}_fixed': str(x_phase.get('fixed_duration')).lower() == 'true'})
                            d_parameter.update({f'{s_celltype}_cycle_phase_transition_rate_{s_index}': float(x_phase.text)})
                            ds_unit.update({f'{s_celltype}_cycle_phase_transition_rate_{s_index}': x_phase.get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cycle><phase_transition_rates> node missing.')
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cycle> node missing.')

                # phenotype death model
                x_death = x_phenotype.find('death')
                if not (x_death is None):
                    # model
                    for x_model in x_death.findall('model'):
                        s_code = str(x_model.get('code'))
                        s_model = ds_death_model[s_code]  # apoptosis oder necrosis
                        d_parameter.update({f'{s_celltype}_death_{s_model}_rate': float(x_model.find('death_rate').text)})
                        # duration
                        try:
                            for x_phase in x_model.findall('duration'):
                                s_index = str(x_phase.get('index'))
                                d_parameter.update({f'{s_celltype}_death_{s_model}_phase_{s_index}_fixed': str(x_phase.get('fixed_duration')).lower() == 'true'})
                                d_parameter.update({f'{s_celltype}_death_{s_model}_phase_{s_index}': float(x_phase.text)})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><phase_durations> node missing.')
                        # rate
                        try:
                            for x_phase in x_model.findall('rate'):
                                s_index = str(x_phase.get('start_index')) + _ + str(x_phase.get('end_index'))
                                d_parameter.update({f'{s_celltype}_death_{s_model}_phase_{s_index}_fixed': str(x_phase.get('fixed_duration')).lower() == 'true'})
                                d_parameter.update({f'{s_celltype}_death_{s_model}_phase_{s_index}': float(x_phase.text)})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><phase_transition_rates> node missing.')
                        # parameters
                        x_death_model_parameter = x_model.find('parameters')
                        if not (x_death_model_parameter is None):
                            # unlysed_fluid_change_rate
                            try:
                                d_parameter.update({f'{s_celltype}_death_{s_model}_unlysed_fluid_change_rate': float(x_death_model_parameter.find('unlysed_fluid_change_rate').text)})
                                ds_unit.update({f'{s_celltype}_death_{s_model}_unlysed_fluid_change_rate': x_death_model_parameter.find('unlysed_fluid_change_rate').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><unlysed_fluid_change_rate> node missing.')
                            # lysed_fluid_change_rate
                            try:
                                d_parameter.update({f'{s_celltype}_death_{s_model}_lysed_fluid_change_rate': float(x_death_model_parameter.find('lysed_fluid_change_rate').text)})
                                ds_unit.update({f'{s_celltype}_death_{s_model}_lysed_fluid_change_rate': x_death_model_parameter.find('lysed_fluid_change_rate').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><lysed_fluid_change_rate> node missing.')
                            # cytoplasmic_biomass_change_rate
                            try:
                                d_parameter.update({f'{s_celltype}_death_{s_model}_cytoplasmic_biomass_change_rate': float(x_death_model_parameter.find('cytoplasmic_biomass_change_rate').text)})
                                ds_unit.update({f'{s_celltype}_death_{s_model}_cytoplasmic_biomass_change_rate': x_death_model_parameter.find('cytoplasmic_biomass_change_rate').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><cytoplasmic_biomass_change_rate> node missing.')
                            # nuclear_biomass_change_rate
                            try:
                                d_parameter.update({f'{s_celltype}_death_{s_model}_nuclear_biomass_change_rate': float(x_death_model_parameter.find('nuclear_biomass_change_rate').text)})
                                ds_unit.update({f'{s_celltype}_death_{s_model}_nuclear_biomass_change_rate': x_death_model_parameter.find('nuclear_biomass_change_rate').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><nuclear_biomass_change_rate> node missing.')
                            # calcification_rate
                            try:
                                d_parameter.update({f'{s_celltype}_death_{s_model}_calcification_rate': float(x_death_model_parameter.find('calcification_rate').text)})
                                ds_unit.update({f'{s_celltype}_death_{s_model}_calcification_rate': x_death_model_parameter.find('calcification_rate').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><calcification_rate> node missing.')
                            # relative_rupture_volume
                            try:
                                d_parameter.update({f'{s_celltype}_death_{s_model}_relative_rupture_volume': float(x_death_model_parameter.find('relative_rupture_volume').text)})
                                ds_unit.update({f'{s_celltype}_death_{s_model}_relative_rupture_volume': x_death_model_parameter.find('relative_rupture_volume').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><relative_rupture_volume> node missing.')
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death><model code={s_code} name={s_model}><parameters> node missing.')
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><death> node missing.')

                # phenotype live volume
                x_volume = x_phenotype.find('volume')
                if not (x_volume is None):
                    # total
                    try:
                        d_parameter.update({f'{s_celltype}_volume_total': float(x_volume.find('total').text)})
                        ds_unit.update({f'{s_celltype}_volume_total': x_volume.find('total').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><total> node missing.')
                    # fluid_fraction
                    try:
                        d_parameter.update({f'{s_celltype}_volume_fluid_fraction': float(x_volume.find('fluid_fraction').text)})
                        ds_unit.update({f'{s_celltype}_volume_fluid_fraction': x_volume.find('fluid_fraction').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><fluid_fraction> node missing.')
                    # nuclear
                    try:
                        d_parameter.update({f'{s_celltype}_volume_nuclear': float(x_volume.find('nuclear').text)})
                        ds_unit.update({f'{s_celltype}_volume_nuclear': x_volume.find('nuclear').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><nuclear> node missing.')
                    # fluid_change_rate
                    try:
                        d_parameter.update({f'{s_celltype}_volume_fluid_change_rate': float(x_volume.find('fluid_change_rate').text)})
                        ds_unit.update({f'{s_celltype}_volume_fluid_change_rate': x_volume.find('fluid_change_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><fluid_change_rate> node missing.')
                    # cytoplasmic_biomass_change_rate
                    try:
                        d_parameter.update({f'{s_celltype}_volume_cytoplasmic_biomass_change_rate': float(x_volume.find('cytoplasmic_biomass_change_rate').text)})
                        ds_unit.update({f'{s_celltype}_volume_cytoplasmic_biomass_change_rate': x_volume.find('cytoplasmic_biomass_change_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><cytoplasmic_biomass_change_rate> node missing.')
                    # nuclear_biomass_change_rate
                    try:
                        d_parameter.update({f'{s_celltype}_volume_nuclear_biomass_change_rate': float(x_volume.find('nuclear_biomass_change_rate').text)})
                        ds_unit.update({f'{s_celltype}_volume_nuclear_biomass_change_rate': x_volume.find('nuclear_biomass_change_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><nuclear_biomass_change_rate> node missing.')
                    # calcified_fraction
                    try:
                        d_parameter.update({f'{s_celltype}_volume_calcified_fraction': float(x_volume.find('calcified_fraction').text)})
                        ds_unit.update({f'{s_celltype}_volume_calcified_fraction': x_volume.find('calcified_fraction').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><calcified_fraction> node missing.')
                    # calcification_rate
                    try:
                        d_parameter.update({f'{s_celltype}_volume_calcification_rate': float(x_volume.find('calcification_rate').text)})
                        ds_unit.update({f'{s_celltype}_volume_calcification_rate': x_volume.find('calcification_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><calcification_rate> node missing.')
                    # relative_rupture_volume
                    try:
                        d_parameter.update({f'{s_celltype}_volume_relative_rupture_volume': float(x_volume.find('relative_rupture_volume').text)})
                        ds_unit.update({f'{s_celltype}_volume_relative_rupture_volume': x_volume.find('relative_rupture_volume').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume><relative_rupture_volume> node missing.')
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><volume> node missing.')

                # phenotype live mechanics
                x_mechanics = x_phenotype.find('mechanics')
                if not (x_mechanics is None):
                    # cell_cell_adhesion_strength
                    try:
                        d_parameter.update({f'{s_celltype}_mechanics_cell_cell_adhesion_strength': float(x_mechanics.find('cell_cell_adhesion_strength').text)})
                        ds_unit.update({f'{s_celltype}_mechanics_cell_cell_adhesion_strength': x_mechanics.find('cell_cell_adhesion_strength').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><cell_cell_adhesion_strength> node missing.')
                    # cell_cell_repulsion_strength
                    try:
                        d_parameter.update({f'{s_celltype}_mechanics_cell_cell_repulsion_strength': float(x_mechanics.find('cell_cell_repulsion_strength').text)})
                        ds_unit.update({f'{s_celltype}_mechanics_cell_cell_repulsion_strength': x_mechanics.find('cell_cell_repulsion_strength').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><cell_cell_repulsion_strength> node missing.')
                    # cell_BM_adhesion_strength
                    try:
                        d_parameter.update({f'{s_celltype}_mechanics_cell_BM_adhesion_strength': float(x_mechanics.find('cell_BM_adhesion_strength').text)})
                        ds_unit.update({f'{s_celltype}_mechanics_cell_BM_adhesion_strength': x_mechanics.find('cell_BM_adhesion_strength').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><cell_BM_adhesion_strength> node missing.')
                    # cell_BM_repulsion_strength
                    try:
                        d_parameter.update({f'{s_celltype}_mechanics_cell_BM_repulsion_strength': float(x_mechanics.find('cell_BM_repulsion_strength').text)})
                        ds_unit.update({f'{s_celltype}_mechanics_cell_BM_repulsion_strength': x_mechanics.find('cell_BM_repulsion_strength').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><cell_BM_repulsion_strength> node missing.')
                    # attachment_elastic_constant
                    try:
                        d_parameter.update({f'{s_celltype}_mechanics_elastic_constant': float(x_mechanics.find('attachment_elastic_constant').text)})
                        ds_unit.update({f'{s_celltype}_mechanics_elastic_constant': x_mechanics.find('attachment_elastic_constant').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><attachment_elastic_constant> node missing.')
                    # attachment_rate
                    try:
                        d_parameter.update({f'{s_celltype}_mechanics_attachment_rate': float(x_mechanics.find('attachment_rate').text)})
                        ds_unit.update({f'{s_celltype}_mechanics_attachment_rate': x_mechanics.find('attachment_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><attachment_rate> node missing.')
                    # detachment_rate
                    try:
                        d_parameter.update({f'{s_celltype}_mechanics_detachment_rate': float(x_mechanics.find('detachment_rate').text)})
                        ds_unit.update({f'{s_celltype}_mechanics_detachment_rate': x_mechanics.find('detachment_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><detachment_rate> node missing.')
                    # relative_maximum_adhesion_distance
                    try:
                        d_parameter.update({f'{s_celltype}_mechanics_relative_maximum_adhesion_distance': float(x_mechanics.find('relative_maximum_adhesion_distance').text)})
                        ds_unit.update({f'{s_celltype}_mechanics_relative_maximum_adhesion_distance': x_mechanics.find('relative_maximum_adhesion_distance').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><relative_maximum_adhesion_distance> node missing.')
                    # options
                    x_mechanics_options = x_mechanics.find('options')
                    if not (x_mechanics_options is None):
                        # set_relative_equilibrium_distance
                        try:
                            d_parameter.update({f'{s_celltype}_mechanics_relative_equilibrium_distance_enabled': str(x_mechanics_options.find('set_relative_equilibrium_distance').get('enabled').lower() == 'true')})
                            d_parameter.update({f'{s_celltype}_mechanics_relative_equilibrium_distance': float(x_mechanics_options.find('set_relative_equilibrium_distance').text)})
                            ds_unit.update({f'{s_celltype}_mechanics_relative_equilibrium_distance': x_mechanics_options.find('set_relative_equilibrium_distance').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><options><set_relative_equilibrium_distance> node missing.')
                        # set_absolute_equilibrium_distance
                        try:
                            d_parameter.update({f'{s_celltype}_mechanics_absolute_equilibrium_distance_enabled': str(x_mechanics_options.find('set_absolute_equilibrium_distance').get('enabled').lower() == 'true')})
                            d_parameter.update({f'{s_celltype}_mechanics_absolute_equilibrium_distance': float(x_mechanics_options.find('set_absolute_equilibrium_distance').text)})
                            ds_unit.update({f'{s_celltype}_mechanics_absolute_equilibrium_distance': x_mechanics_options.find('set_absolute_equilibrium_distance').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><options><set_absolute_equilibrium_distance> node missing.')
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><options> node missing.')
                    # cell_adhesion_affinities
                    x_mechanics_affinities = x_mechanics.find('cell_adhesion_affinities')
                    if not (x_mechanics_affinities is None):
                        for x_affinity in x_mechanics_affinities.findall('cell_adhesion_affinity'):
                            s_name = x_affinity.get('name')
                            d_parameter.update({f'{s_celltype}_mechanics_relative_cell_adhesion_affinity_{s_name}': float(x_affinity.text)})
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics><cell_adhesion_affinities> node missing.')
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><mechanics> node missing.')

                # phenotype live motility
                x_motility = x_phenotype.find('motility')
                if not (x_mechanics is None):
                    # speed
                    try:
                        d_parameter.update({f'{s_celltype}_motility_speed': float(x_motility.find('speed').text)})
                        ds_unit.update({f'{s_celltype}_motility_speed': x_motility.find('speed').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><motility><speed> node missing.')
                    # persistence_time
                    try:
                        d_parameter.update({f'{s_celltype}_motility_persistence_time': float(x_motility.find('persistence_time').text)})
                        ds_unit.update({f'{s_celltype}_motility_persistence_time': x_motility.find('persistence_time').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><motility><persistence_time> node missing.')
                    # migration_bias
                    try:
                        d_parameter.update({f'{s_celltype}_motility_migration_bias': float(x_motility.find('migration_bias').text)})
                        ds_unit.update({f'{s_celltype}_motility_migration_bias': x_motility.find('migration_bias').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><motility><_motility_migration_bias> node missing.')
                    # option
                    x_motility_options = x_motility.find('options')
                    if not (x_motility_options is None):
                        d_parameter.update({f'{s_celltype}_motility_enabled': str(x_motility_options.find('enabled').text).lower() == 'true'})
                        d_parameter.update({f'{s_celltype}_motility_use_2D': str(x_motility_options.find('use_2D').text).lower() == 'true'})
                        # optionn chemotaxis
                        x_motility_chemotaxis = x_motility.find('options').find('chemotaxis')
                        if not (x_motility_chemotaxis):
                            d_parameter.update({f'{s_celltype}_motility_chemotaxis_enabled': str(x_motility_chemotaxis.find('enabled').text).lower() == 'true'})
                            d_parameter.update({f'{s_celltype}_motility_chemotaxis_substrate': str(x_motility_chemotaxis.find('substrate').text)})
                            d_parameter.update({f'{s_celltype}_motility_chemotaxis_direction': float(x_motility_chemotaxis.find('direction').text)})
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><motility><options><chemotaxis> node missing.')
                        # options advanced_chemotaxis
                        x_motility_advancedchemotaxis = x_motility.find('options').find('advanced_chemotaxis')
                        if not (x_motility_advancedchemotaxis is None):
                            d_parameter.update({f'{s_celltype}_motility_advanced_chemotaxis_enabled': str(x_motility_advancedchemotaxis.find('enabled').text).lower() == 'true'})
                            d_parameter.update({f'{s_celltype}_motility_advanced_chemotaxis_normalize_each_gradient': str(x_motility_advancedchemotaxis.find('normalize_each_gradient').text).lower() == 'true'})
                            for x_chemotactic in x_motility_advancedchemotaxis.find('chemotactic_sensitivities').findall('chemotactic_sensitivity'):
                                s_subs = str(x_chemotactic.get('substrate'))
                                d_parameter.update({f'{s_celltype}_motility_advanced_chemotaxis_chemotactic_sensitivity_{s_subs}': float(x_chemotactic.text)})
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><motility><options><advanced_chemotaxis> node missing.')
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><motility><options> node missing.')
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><motility> node missing.')

                # phenotype live secretion
                x_secretion = x_phenotype.find('secretion')
                if not (x_secretion is None):
                    for x_substrate in x_secretion.findall('substrate'):
                        s_subs = str(x_substrate.get('name'))
                        # secretion_rate
                        try:
                            d_parameter.update({f'{s_celltype}_{s_subs}_secretion_rate': float(x_substrate.find('secretion_rate').text)})
                            ds_unit.update({f'{s_celltype}_{s_subs}_secretion_rate': x_substrate.find('secretion_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><secretion><substrate><secretion_rate> node missing.')
                        # secretion_target
                        try:
                            d_parameter.update({f'{s_celltype}_{s_subs}_secretion_target': float(x_substrate.find('secretion_target').text)})
                            ds_unit.update({f'{s_celltype}_{s_subs}_secretion_target': x_substrate.find('secretion_target').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><secretion><substrate><secretion_target> node missing.')
                        # uptake_rate
                        try:
                            d_parameter.update({f'{s_celltype}_{s_subs}_uptake_rate': float(x_substrate.find('uptake_rate').text)})
                            ds_unit.update({f'{s_celltype}_{s_subs}_uptake_rate': x_substrate.find('uptake_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><secretion><substrate><uptake_rate> node missing.')
                        # net_export_rate
                        try:
                            d_parameter.update({f'{s_celltype}_{s_subs}_net_export_rate': float(x_substrate.find('net_export_rate').text)})
                            ds_unit.update({f'{s_celltype}_{s_subs}_net_export_rate': x_substrate.find('net_export_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><secretion><substrate><net_export_rate> node missing.')
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><secretion> node missing.')

                # phenotype live cell_interactions
                x_interaction = x_phenotype.find('cell_interactions')
                if not (x_interaction is None):
                    # dead_phagocytosis_rate
                    try:
                        d_parameter.update({f'{s_celltype}_cell_interaction_dead_phagocytosis_rate': float(x_interaction.find('dead_phagocytosis_rate').text)})
                        ds_unit.update({f'{s_celltype}_cell_interaction_dead_phagocytosis_rate': x_interaction.find('dead_phagocytosis_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cell_interactions><dead_phagocytosis_rate> node missing.')
                    # damage_rate
                    try:
                        d_parameter.update({f'{s_celltype}_cell_interaction_damage_rate': float(x_interaction.find('damage_rate').text)})
                        ds_unit.update({f'{s_celltype}_cell_interaction_damage_rate': x_interaction.find('damage_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cell_interactions><damage_rate> node missing.')
                    # live_phagocytosis_rates
                    x_interaction_phagocytosis = x_interaction.find('live_phagocytosis_rates')
                    if not (x_interaction_phagocytosis is None):
                        for x_rate in x_interaction_phagocytosis:
                            s_name = x_rate.get('name')
                            d_parameter.update({f'{s_celltype}_cell_interaction_{s_name}_live_phagocytosis_rate': float(x_rate.text)})
                            ds_unit.update({f'{s_celltype}_cell_interaction_{s_name}_live_phagocytosis_rate': x_rate.get('units')})
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cell_interactions><live_phagocytosis_rates> node missing.')
                    # attack_rates
                    x_interaction_attack = x_interaction.find('attack_rates')
                    if not (x_interaction_attack is None):
                        for x_rate in x_interaction_attack:
                            s_name = x_rate.get('name')
                            d_parameter.update({f'{s_celltype}_cell_interaction_{s_name}_attack_rate': float(x_rate.text)})
                            ds_unit.update({f'{s_celltype}_cell_interaction_{s_name}_attack_rate': x_rate.get('units')})
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cell_interactions><attack_rates> node missing.')
                    # fusion_rates
                    x_interaction_fusion =  x_interaction.find('fusion_rates')
                    if not (x_interaction_fusion is None):
                        for x_rate in x_interaction_fusion:
                            s_name = x_rate.get('name')
                            d_parameter.update({f'{s_celltype}_cell_interaction_{s_name}_fusion_rate': float(x_rate.text)})
                            ds_unit.update({f'{s_celltype}_cell_interaction_{s_name}_fusion_rate': x_rate.get('units')})
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cell_interactions><fusion_rates> node missing.')

                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cell_interactions> node missing.')

                # phenotype live cell_transformation
                x_transformation = x_phenotype.find('cell_transformations')
                if not (x_transformation is None):
                    for x_rate  in x_transformation.find('transformation_rates').findall('transformation_rate'):
                        s_name = str(x_rate.get('name'))
                        d_parameter.update({f'{s_celltype}_{s_name}_transformation_rate': float(x_rate.text)})
                        ds_unit.update({f'{s_celltype}_{s_name}_transformation_rate': x_rate.get('units')})
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype><cell_transformations> node missing.')

            # phenotype
            else:
                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><phenotype> node missing.')

            # search for custom data
            try:
                for x_element in x_celltype.find('custom_data').iter():
                    if (x_element.tag != 'custom_data'):
                        try:
                            d_parameter.update({x_element.tag : self.custom_type[x_element.tag](x_element.text)})
                        except KeyError:
                            d_parameter.update({x_element.tag: float(x_element.text)})
                            es_customdata.add(x_element.tag)
            except AttributeError:
                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{i_id}"><custom_data> node missing.')

        # if custom data was found
        if (len(es_customdata) > 0):
            print(f'Warning @ pyMCDS._read_setting_xml : cell_definition custom_data without variable type setting detected. {sorted(es_customdata)}')

        ### find the user_parameters node ###
        for x_element in x_root.find('user_parameters').iter():
            if (x_element.tag != 'user_parameters'):
                s_type = x_element.get('type')
                if s_type in {'double', 'real', 'float'}:
                    o_type = float
                elif s_type in {'int', 'integer'}:
                    o_type = int
                elif s_type in {'str', 'string'}:
                    o_type = str
                elif s_type in {'bool','boole'}:
                    print(f'Warning @  pyMCDS._read_setting.xml: user_parameter bool variable type detected. value {x_element.tag} might be ambigiuos. translate into string.')
                    o_type = str
                else:
                    print(f'Warning @  pyMCDS._read_setting.xml: user_parameter unknown variable type {s_type} detected. will translate into string.')
                    o_type = str
                d_parameter.update({x_element.tag: o_type(x_element.text)})

        ### find the cell_rules node ###
        x_rule = x_root.find('cell_rules')
        if not (x_rule is None):
            for x_ruleset in x_rule.find('rulesets'):
                s_pathfile = s_path + '/' + x_ruleset.find('filename').text
                try:
                    df_rule = pd.read_csv(s_pathfile, sep=',')
                    df_rule.columns = ['cell_type','signal','direction','behavoir','saturation_value','base_value','half_max','hill_power','apply_to_dead']
                    if (df_ruleset is None):
                        df_ruleset = df_rule
                    else:
                        df_ruleset = pd.concate([df_ruleset, df_rule], axis=0)
                except pd._libs.parsers.EmptyDataError:
                    print(f'Warning @ pyMCDS._read_setting_xml : {s_pathfile} is empty.')
        else:
            print(f'Warning @ pyMCDS._read_setting_xml : <cell_rules> node missing.')

        ### output ###
        # generate setting dictionary and store extraction from PhysiCell_settings.xml
        d_mcds['metadata'] = {}
        d_mcds['metadata']['substrate'] = ds_substrate
        d_mcds['metadata']['cell_type'] = ds_celltype
        d_mcds['setting'] = {}
        d_mcds['setting']['parameters'] = d_parameter
        d_mcds['setting']['units'] = ds_unit
        d_mcds['setting']['rules'] = df_ruleset
        print('done!')
        return d_mcds

