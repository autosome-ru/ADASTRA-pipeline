import os
import pathlib
component_dict = [
    'alignments_path',
    'scripts_path',
    'configs_path',
    'parallel_parameters_path',
    'badmaps_path',
    'results_path',
    'badmaps_dict_path',
    'master_list_path',
    'tf_dict_path',
    'cl_dict_path',
    'reference_path',
    'FA',
    'repeats_path',
    'genome_path',
    'intervals_path',
    'dbsnp_vcf_path',
]


def remove_around_punctuation(string, with_quotes=False):
    string = string.strip()
    if with_quotes:
        if string and (string.startswith("\"") or string.startswith("'")):
            string = string[1:]
            if string.endswith("\"") or string.endswith("'"):
                string = string[:-1]
            else:
                raise AssertionError('Quotes are not closed', string)
        else:
            raise AssertionError('Value is without quotation', string)
    return string.strip()


def parse_line(line):
    line = line.split('=')
    if len(line) != 2:
        raise AssertionError('Wrong format cfg file!', line)
    return remove_around_punctuation(line[0]), remove_around_punctuation(line[1], True)


def construct_line(component_name, component_value):
    return "{}='{}'\n".format(component_name, component_value)


def pack_line(config_dict, component_name):
    if component_name == 'reference_path':
        return construct_line(component_name, os.path.join(config_dict['scripts_path'], 'Configs', 'reference'))
    elif component_name == 'alignments_path':
        return construct_line(component_name, config_dict[component_name])
    elif component_name == 'results_path':
        return construct_line(component_name, config_dict[component_name])
    elif component_name == 'badmaps_path':
        return construct_line(component_name, config_dict[component_name])
    elif component_name == 'FA':
        return construct_line(component_name, os.path.join(config_dict['scripts_path'], 'Configs', 'genome-norm.fasta'))
    elif component_name == 'badmaps_dict_path':
        return construct_line(component_name,
                              os.path.join(config_dict['scripts_path'], 'Configs', 'badmaps_dict.json'))
    elif component_name == 'cl_dict_path':
        return construct_line(component_name,
                              os.path.join(config_dict['scripts_path'], 'Configs', 'cl_dict.json'))
    elif component_name == 'tf_dict_path':
        return construct_line(component_name,
                              os.path.join(config_dict['scripts_path'], 'Configs', 'tf_dict.json'))
    elif component_name == 'master_list_path':
        return construct_line(component_name, config_dict[component_name])
    elif component_name == 'dbsnp_vcf_path':
        return construct_line(component_name, config_dict[component_name])
    elif component_name == 'scripts_path':
        return construct_line(component_name, config_dict['scripts_path'])
    elif component_name == 'configs_path':
        return construct_line(component_name, os.path.join(config_dict['scripts_path'], 'Configs'))
    elif component_name == 'parallel_parameters_path':
        return construct_line(component_name, os.path.join(config_dict['scripts_path'], 'Configs',
                                                           'parallel_configs'))
    elif component_name == 'repeats_path':
        return construct_line(component_name, config_dict[component_name])
    elif component_name == 'intervals_path':
        return construct_line(component_name, config_dict[component_name])
    elif component_name == 'genome_path':
        return construct_line(component_name, config_dict[component_name])
    raise AssertionError(component_name, 'is not in valid arguments, check Config.cfg file')


def read_cfg_file(cfg_file):
    config_dict = {}
    iteration = ''

    with open(cfg_file) as cfg_buffer:
        for line in cfg_buffer:
            if line.strip():
                if not line.startswith('#'):
                    if iteration == 'soft' or iteration == 'parameters':
                        continue
                    config, value = parse_line(line)
                    config_dict[config] = value
                else:
                    iteration = line.strip()[1:]
    config_dict['scripts_path'] = os.path.join(pathlib.Path(__file__).parent.absolute(), 'scripts')
    with open(os.path.join(config_dict['scripts_path'],
                           'HELPERS', 'paths_for_components.py'), 'w') as out:
        for component_name in component_dict:
            out.write(pack_line(config_dict, component_name))


if __name__ == '__main__':
    read_cfg_file(os.path.join(pathlib.Path(__file__).parent.absolute(), 'CONFIG.cfg'))
