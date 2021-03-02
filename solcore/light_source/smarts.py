import numpy as np
import os
import shutil
import subprocess
import tempfile
import platform
import solcore
from datetime import datetime


class SmartsSolverError(Exception):
    pass


smarts = solcore.config.smarts()


system = platform.system()
if system == 'Windows':
    extension = '.exe'
elif system == 'Linux':
    extension = ''
else:
    extension = '.command'

executable = os.path.join(smarts, "smarts295" + extension)

error_msg = """ERROR: SMARTS location not correctly configured or SMARTS executable not working.

SMARTS can be obtained free of charge from the NREL webpage:

http://www.nrel.gov/rredc/smarts/

You might need to re-compile the Fortran code as current 64 bit CPUs are not supported by the shipped binaries.
"""


def skipper(fname):
    with open(fname) as fin:
        no_comments = (line for line in fin if "   " in line)
        next(no_comments, None)  # skip header
        next(no_comments, None)  # skip header
        for row in no_comments:
            yield row


def calculate_spectrum_smarts(smarts_file_contents=None, filename='smarts295', target_directory=None):
    """

    :param smarts_file_contents:
    :param filename:
    :param target_directory:
    :return:
    """
    if not smarts:
        raise SmartsSolverError(f"Smarts installation not found in {smarts}")

    if smarts_file_contents is None:
        smarts_file_contents = build_smarts_file(**get_default_smarts_object())
    else:
        smarts_file_contents = build_smarts_file(**smarts_file_contents)

    if target_directory is not None:
        target_directory = target_directory.rstrip('\/')
        assert os.access(target_directory, os.W_OK), 'ERROR: Target folder for smarts output does not exists ' \
                                                     'or is not writable.'

    data = []
    # We use a temp directory to store temporary the data.
    # If needed, we'll copy it later to the target directory.
    with tempfile.TemporaryDirectory(prefix="tmp", suffix="_sc3SMARTS") as working_directory:

        ext_file = os.path.join(working_directory, "{}.ext.txt".format(filename))
        inp_file = os.path.join(working_directory, "{}.inp.txt".format(filename))
        out_file = os.path.join(working_directory, "{}.out.txt".format(filename))
        scn_file = os.path.join(working_directory, "{}.scn.txt".format(filename))

        # Save the data to the input file
        with open(inp_file, "w") as f:
            f.write(smarts_file_contents)

        # Start the process
        try:
            this_process = subprocess.Popen((executable,), stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                            stderr=subprocess.PIPE, cwd=smarts)
            # We need to tell smarts where to find the input data
            output, error = this_process.communicate(
                input=bytes('N\n"{0}"\n{1}\nY\n'.format(working_directory, filename), "ASCII"))

            error = error.decode('utf8')

            if len(error) > 0 and 'Note' not in error:
                raise RuntimeError(error)

            data = []
            if os.path.exists(scn_file):
                data = np.loadtxt(skipper(scn_file), unpack=True)
                if target_directory is not None:
                    shutil.copy2(scn_file, scn_file.replace(working_directory, target_directory))
            if os.path.exists(inp_file):
                if target_directory is not None:
                    shutil.copy2(inp_file, inp_file.replace(working_directory, target_directory))
            if os.path.exists(ext_file):
                if target_directory is not None:
                    shutil.copy2(ext_file, ext_file.replace(working_directory, target_directory))
            if os.path.exists(out_file):
                if target_directory is not None:
                    shutil.copy2(out_file, out_file.replace(working_directory, target_directory))

        except RuntimeError as err:
            print('ERROR in SMARTS: {}'.format(err))

        except ValueError as err:
            print('ERROR in SMARTS: {}'.format(err))

        except Exception as err:
            print('ERROR in SMARTS: {}'.format(error_msg))
            print('ERROR in SMARTS: {}'.format(err))

    if len(data) == 0:
        raise ValueError('ERROR in SMARTS: Output file is empty. Likely error in the input parameters.')

    return data


def build_smarts_file(**kwargs):
    try:
        smarts_file_contents = ["'{COMNT}' !Card 1"]
        smarts_file_contents.append("{ISPR} !Card 2")
        smarts_file_contents.append([
                                        "{SPR} !Card 2a",
                                        "{SPR} {ALTIT} {HEIGHT} !Card 2a",
                                        "{LATIT} {ALTIT} {HEIGHT} !Card 2a"
                                    ][kwargs["ISPR"]])

        smarts_file_contents.append("{IATMOS} !Card 3")
        smarts_file_contents.append([
                                        "{TAIR} {RH} '{SEASON}' {TDAY} !Card 3a",
                                        "'{ATMOS}' !Card 3a"
                                    ][kwargs["IATMOS"]])

        smarts_file_contents.append("{IH2O} !Card 4")
        if kwargs["IH2O"] == 0:
            smarts_file_contents.append("{W} !Card 4a")

        smarts_file_contents.append("{IO3} !Card 5")
        if kwargs["IO3"] == 0:
            smarts_file_contents.append("{IALT} {AbO3} !Card 5a")

        smarts_file_contents.append("{IGAS} !Card 6")
        if kwargs["IGAS"] == 0:
            smarts_file_contents.append("{ILOAD} !Card 6a")
            if kwargs["ILOAD"] == 0:
                smarts_file_contents.append(
                    "{ApCH2O} {ApCH4} {ApCO} {ApHNO2} {ApHNO3} {ApNO} {ApNO2} {ApNO3} {ApO3} {ApSO2} !Card 6b"
                )

        smarts_file_contents.append("{qCO2} !Card 7")
        smarts_file_contents.append("{ISPCTR} !Card 7a")

        smarts_file_contents.append("'{AEROS}' !Card 8")
        if kwargs["AEROS"] == "USER":
            smarts_file_contents.append("{ALPHA1} {ALPHA2} {OMEGL} {GG} !Card 8a")

        smarts_file_contents.append("{ITURB} !Card 9")
        smarts_file_contents.append([
                                        "{TAU5} !Card 9a",
                                        "{BETA} !Card 9a",
                                        "{BCHUEP} !Card 9a",
                                        "{RANGE} !Card 9a",
                                        "{VISI} !Card 9a",
                                        "{TAU550} !Card 9a",
                                    ][kwargs["ITURB"]])

        smarts_file_contents.append("{IALBDX} !Card 10")
        if kwargs["IALBDX"] == -1:
            smarts_file_contents.append("{RHOX} !Card 10a")
        smarts_file_contents.append("{ITILT} !Card 10b")
        if kwargs["ITILT"] == 1:
            smarts_file_contents.append("{IALBDG} {TILT} {WAZIM} !Card 10c")
            if kwargs["IALBDG"] == 1:
                smarts_file_contents.append("{RHOG} !Card 10d")

        smarts_file_contents.append("{WLMN} {WLMX} {SUNCOR} {SOLARC} !Card 11")

        smarts_file_contents.append("{IPRT} !Card 12")
        if kwargs["IPRT"] >= 1:
            smarts_file_contents.append("{WPMN} {WPMX} {INTVL} !Card 12a")
            if kwargs["IPRT"] == 2 or kwargs["IPRT"] == 3:
                smarts_file_contents.append("{IOTOT} !Card 12b")
                smarts_file_contents.append("{IOUT} !Card 12c")

        smarts_file_contents.append("{ICIRC} !Card 13")
        if kwargs["ICIRC"] == 1:
            smarts_file_contents.append("{SLOPE} {APERT} {LIMIT} !Card 13a")

        smarts_file_contents.append("{ISCAN} !Card 14")
        if kwargs["ISCAN"] == 1:
            smarts_file_contents.append("{IFILT} {WV1} {WV2} {STEP} {FWHM} !Card 14a")

        smarts_file_contents.append("{ILLUM} !Card 15")

        smarts_file_contents.append("{IUV} !Card 16")

        smarts_file_contents.append("{IMASS} !Card 17")
        smarts_file_contents.append([
                                        "{ZENIT} {AZIM} !Card 17a",
                                        "{ELEV} {AZIM} !Card 17a",
                                        "{AMASS} !Card 17a",
                                        "{YEAR} {MONTH} {DAY} {HOUR} {LATIT} {LONGIT} {ZONE} !Card 17a",
                                        "{MONTH} {LATIT} {DSTEP} !Card 17a",
                                    ][kwargs["IMASS"]])

        smarts_file_schema = "\n".join(smarts_file_contents)
        # print (smarts_file_schema)

        smarts_file_complete = smarts_file_schema.format(**kwargs) + "\n"
    except KeyError as err:
        print("The SMARTS options you have selected require additional data, variables are undefined.")
        print(err)
        raise

    return smarts_file_complete


def get_default_smarts_object():
    """  Provides a dictionary with most of the parameters (values of the CARDS) needed by SMARTS to calculate
     a solar spectrum. It can be used as a template to customize for user defined conditions.

    :return:
    """

    # CONSTANT PARAMETERS   ---------------------------------
    # location and general time info
    latitude = 40.4966  # 'LATIT', deg, latitude
    longitude = -3.4620  # 'LONGIT, deg, longitude
    preasure_model = 1  # 'ISPR', surface preasure model set to real data + altitude correction
    altitude = 0.625  # 'ALTIT', km, altitude above sea level
    altura = 0.0  # 'HEIGHT', m, altude over the ground level
    time_zone = 0  # 'ZONE', time zone
    season = 'SUMMER'  # 'SEASON', season of the year. Can be summer or winter only
    albedo = 9  # 'IALBDX', 'IALBDG', Albedo model = 9, Dry clay soil
    solar_position_mode = 3  # 'IMASS', use the location and time to calculate the position of the sun

    # atmospheric conditions
    atmospheric_data = 1  # 'IATMOS', allows to input the correct atmospheric data
    atmosphere_model = 'USSA'  # 'ATMOS', 'US standard spectrum', general atmospheric model, setting all parameters not given as input
    precipitable_water = 0  # 'IH2O', water vapor data given as input
    ozone = 1  # 'IO3', ozone abundance calculated from atmospheric data
    gas_contents = 0  # 'IGAS', gas abundances set to a particular scenario
    gas_scenairo = 2  # 'ILOAD', light pollution scenario
    CO2content = 370  # 'qCO2', CO2 content in ppmv
    extraterrestial = 0  # 'ISPCTR', extraterrestial spectrum
    aerosol = 'S&F_RURAL'  # 'AEROS', aerosol model
    turbidity = 0  # 'ITURB', select turbidity model
    tau500_param = 0.085  # 'TAU5'

    # output info
    print_info = 2  # 'IPRT', print output in spreadsheet format
    total_variables = 1  # 'IOTOT', total number of variables to print
    which_variables = '4'  # 'IOUT', code of the output variables. Set to all tilted irradiances + experimental with FoV
    wavelenght_min = 280  # 'WLMN', 'WPMN', nm
    wavelenght_max = 4004  # 'WLMX', 'WPMX', nm
    wavelenght_step = 0.5  # 'INTVL', nm
    convolute = 1  # 'ISCAN'
    conv_function = 1  # 'IFILT', Gausian
    conv_wl_min = 300  # 'WV1', nm
    conv_wl_max = 3990  # 'WV2', nm
    conv_FWHM = 4  # 'FWHM', nm
    conv_step = 2  # 'STEP', nm
    tilt = 1  # 'ITILT', enable calculations for a tilt surface
    altitude_tilt = -999  # 'TILT', -999 means 'track the Sun'
    azimuth_tilt = -999  # 'WAZIM', -999 means 'track the Sun'
    solar_constant = 1367  # 'SOLARC, W/m2, solar constant
    sun_correction = 1  # 'SUNCOR', distance to the Sun correction faction. It is not used

    # Colimator information
    circumsolar = 1  # 'ICIRC', if circumsolar radiation is to be calculated
    slope = 1  # 'SLOPE', configuration of the colimator
    aperture = 2.5  # 'APERT'
    limit = 4  # 'LIMIT'

    # Others
    illuminance = 0  # 'ILLUM'
    UVcalc = 0  # 'IUV'
    air_mass = 3  # 'IMASS', set to calculate the air mass from the time and location

    # VARIABLES     ---------------------------------
    comment = 'Test'  # 'COMNT', single line with a comment on the run. Max 64 characters, no spaces but underscores
    P = 940  # 'SPR', mb, surface preasure
    T_air = 17  # 'TAIR', deg, temperature of the air
    T_day = 25  # 'TDAY', deg, average daily temperature
    humid = 30  # 'RH', %, relative humidity
    water_vapour = 1  # 'W', cm, precipitable water
    targetTime = datetime(2015, 5, 19, 12,
                          30)  # has to be split in 'YEAR', 'MONTH', 'DAY', 'HOUR', the later in Local Standard Time

    smarts_input = {
        'LATIT': latitude,
        'LONGIT': longitude,
        'ISPR': preasure_model,
        'ALTIT': altitude,
        'HEIGHT': altura,
        'ZONE': time_zone,
        'SEASON': season,
        'IALBDX': albedo,
        'IALBDG': albedo,
        'IMASS': solar_position_mode,
        'IATMOS': atmospheric_data,
        'ATMOS': atmosphere_model,
        'IH2O': precipitable_water,
        'IO3': ozone,
        'IGAS': gas_contents,
        'ILOAD': gas_scenairo,
        'qCO2': CO2content,
        'ISPCTR': extraterrestial,
        'AEROS': aerosol,
        'ITURB': turbidity,
        'TAU5': tau500_param,
        'IPRT': print_info,
        'IOTOT': total_variables,
        'IOUT': which_variables,
        'WLMN': wavelenght_min,
        'WPMN': wavelenght_min,
        'WLMX': wavelenght_max,
        'WPMX': wavelenght_max,
        'INTVL': wavelenght_step,
        'ISCAN': convolute,
        'IFILT': conv_function,
        'WV1': conv_wl_min,
        'WV2': conv_wl_max,
        'FWHM': conv_FWHM,
        'STEP': conv_step,
        'ITILT': tilt,
        'TILT': altitude_tilt,
        'WAZIM': azimuth_tilt,
        'SOLARC': solar_constant,
        'SUNCOR': sun_correction,
        'ICIRC': circumsolar,
        'SLOPE': slope,
        'APERT': aperture,
        'LIMIT': limit,
        'ILLUM': illuminance,
        'IUV': UVcalc,
        'COMNT': comment,
        'SPR': P,
        'TAIR': T_air,
        'TDAY': T_day,
        'RH': humid,
        'W': water_vapour,
        'YEAR': targetTime.year,
        'MONTH': targetTime.month,
        'DAY': targetTime.day,
        'HOUR': targetTime.hour + (targetTime.minute + targetTime.second / 60) / 60
    }

    return smarts_input


if __name__ == "__main__":
    data = calculate_spectrum_smarts(target_directory='/Users/diego/Downloads')

    import matplotlib.pyplot as plt

    plt.plot(data[0], data[1])
    plt.plot(data[0], data[2])
    plt.plot(data[0], data[3])
    plt.plot(data[0], data[4])
    plt.show()
