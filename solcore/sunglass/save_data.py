import numpy as np
import pandas as pd
from collections import OrderedDict


def save_iv_data(root_filename, output, light_iv=False):
    """

    :return:
    """
    suffix = 'light' if light_iv else 'dark'
    filename = '{}_{}_{}.csv'.format(root_filename, 'IV', suffix)

    out = OrderedDict({'J': -output.iv['IV'][1], 'V': output.iv['IV'][0]})
    v = output.iv['IV'][0]

    # We get the IV of the normal junctions when included in the MJ structure
    for j in range(len(output.iv['junction IV'])):
        new_j = output.junction_indices[j]
        out['V_MJ{}'.format(new_j)] = output.iv['junction IV'][j][0]

    # And the series resistance
    out['V_Rs'] = output.iv['Rseries IV'][0]

    # We get the IV of the normal junctions if they were isolated
    for j in output.junction_indices:
        out['J_{}'.format(j)] = output[j].iv(v)

    # And the tunnel junctions
    for j in output.tunnel_indices:
        out['J_TJ{}'.format(j)] = output[j].iv(v)

    # Save the data as CSV file
    df = pd.DataFrame(out)
    df.to_csv(filename)
    return df, filename


def save_qe_data(root_filename, output):
    """

    :return:
    """

    qe_data = OrderedDict({'WL': output.wavelength,
                           'R': output.reflected})

    for j in range(len(output)):

        # All layers and junctions have defined a fraction of absorbed light, except some TJ and 2D kind of
        # junctions. We try to get that, first.
        try:
            new_key = 'A_{}'.format(j)
            z = np.arange(0, output[j].width, 1e-11)
            absorbed_per_wl = np.trapz(output[j].absorbed(z), z, axis=0)
            qe_data[new_key] = absorbed_per_wl

        except KeyError:
            print('WARNING: Element {} has no absorption defined.'.format(j))

        # We try to get the QE. If we have a single layer or TJ, this will fail as they don't have QE
        try:
            for k in output[j].qe.keys():
                if k != 'WL':
                    new_key = '{}_{}'.format(k, j)
                    qe_data[new_key] = output[j].qe[k]
        except AttributeError:
            print('WARNING: No "qe" calculated for element {}.'.format(j))

    filename = '{}_{}.csv'.format(root_filename, 'QE')
    df = pd.DataFrame.from_dict(qe_data)
    df.to_csv(filename)
    return df, filename


def save_bandstructure_data(root_filename, output, task):
    """

    :return:
    """
    bs_data = OrderedDict({})
    tsk = '{}_data'.format(task)

    for j in output.junction_indices:

        # Only junctions have information on bandstructure, and not all the models, so we take that into account.
        try:
            data = output[j].__dict__[tsk]['Bandstructure']
            x = data['x']

            if task == 'equilibrium':
                # The properties and the bandstructure are mapped at different positions, so we interpolate
                for p in output[j].__dict__[tsk]['Properties'].keys():
                    old_x = output[j].__dict__[tsk]['Properties']['x']
                    if p != 'x':
                        data[p] = np.interp(x, old_x, output[j].__dict__[tsk]['Properties'][p])
        except KeyError as err:
            print('WARNING: Junction {} has no bandstructure information.'.format(j))
            continue

        # So far, so good: we continue getting the data
        data['x'] += output[j].offset

        new_data = OrderedDict({})
        # We change the key names to include the junction number
        for key in data:
            new_key = '{}_{}'.format(key, j)
            new_data[new_key] = data[key]

        bs_data.update(new_data)

    # If some data has been added to the dictionary, we save that data.
    if len(bs_data.keys()) > 0:
        filename = '{}_{}.csv'.format(root_filename, task)
        df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in bs_data.items()]))
        df.to_csv(filename, index=False)
        return df, filename
