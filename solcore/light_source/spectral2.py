import numpy
from numpy import *
from datetime import datetime
from solcore import spectral_conversion_nm_ev
from solcore.science_tracker import science_reference
import os

this_dir = os.path.split(__file__)[0]


def equation_of_time(day_angle):
    value = 229.18 * (0.000075 + 0.001868 * cos(day_angle) - 0.032077 * sin(day_angle) -
                      0.014615 * cos(2 * day_angle) - 0.040849 * sin(2 * day_angle))
    return value


def loadUtilitySpectra():
    global am_zero_wavelength, am_zero_irradiance, waterspectra, ozonespectra, uniformgasspectra

    spctra = numpy.loadtxt(os.path.join(this_dir, "SPCTRAL_si_units.txt"), unpack=True)
    (am_zero_wavelength, am_zero_irradiance, waterspectra, ozonespectra, uniformgasspectra) \
        = spctra
    waterspectra = waterspectra / 100
    ozonespectra = ozonespectra / 100


loadUtilitySpectra()


def get_default_spectral2_object():
    """ Returns a default state object for use with calculate_spectrum. Parameters returned are for:
        - Coventry (52.39°N,-1.56°E)
        - for 12:30 pm on the 30th of June 2011 (that's in the future!) (not any more..)
        - rural atmospheric conditions
        - 30% humidity
        - 1.42mm precipitated water
        - 1030 hPa atmospheric pressure
        - 0.34mm of Ozone
        - turbidity of 15%.
    """
    defaults = {}

    defaults["latitude"] = 52.39 / 180 * numpy.pi  # Coventry
    defaults["longitude"] = -1.56 / 180 * numpy.pi
    defaults["dateAndTime"] = datetime(2011, 6, 30, 12, 30)  # 12:30pm, 30th of June 2011
    defaults["aod_model"] = "rural"
    defaults["pressure"] = 103000.0  # si("1030 hPa")  #si ("1 atm")
    defaults["humidity"] = 0.3  # si("30 %")
    defaults["precipwater"] = 0.00142  # si("1.42 mm")    # unit?
    defaults["ozone"] = 0.00034  # si("0.34 mm")          # unit?
    defaults["turbidity"] = 0.15  # unit?

    defaults["display units"] = {
        "pressure": ("hPa", 0),
        "latitude": (u"degrees", 3),
        "longitude": (u"degrees", 3),
        "precipwater": ("cm", 3),
        "ozone": ("cm", 3),
        "humidity": ("%", 1),
    }

    return defaults


def calculate_spectrum_spectral2(stateObject=None, suppress_nan=True, power_density_in_nm=False):
    """ Calculates a solar spectrum using the SPECTRAL2 irradiance model developed by the NRL:

    "http://rredc.nrel.gov/solar/models/spectral/SPCTRAL2/"
    
    :param stateObject:
    :param suppress_nan:
    :return:
    """
    science_reference("spectral2 irradiance model", "http://rredc.nrel.gov/solar/models/spectral/SPCTRAL2/")

    if stateObject == None:
        stateObject = get_default_spectral2_object()

    latitude = stateObject["latitude"]
    longitude = stateObject["longitude"]
    dateAndTime = stateObject["dateAndTime"]
    aod_model = stateObject["aod_model"]
    pressure = stateObject["pressure"]
    humidity = stateObject["humidity"]
    ozone = stateObject["ozone"]
    turbidity = stateObject["turbidity"]
    precipwater = stateObject["precipwater"]

    assert aod_model in "rural urban maritime tropospheric".split(), "aod_model must be rural, urban, maritime, or tropospheric"

    time_since_new_years = dateAndTime - datetime(dateAndTime.year, 1, 1, 0, 0)  # 1/1/year, 0:00 am.
    time_since_midnight = dateAndTime - datetime(dateAndTime.year, dateAndTime.month, dateAndTime.day, 0, 0)
    # hours_since_midnight = time_since_midnight.seconds / convert(time_since_midnight.seconds, "s", "h")
    hours_since_midnight = time_since_midnight.seconds / 3600
    day_number = time_since_new_years.days

    # converting back go degrees for longitude rounding
    longitude_degrees = longitude / numpy.pi * 180  # convert(longitude,"radians",u"degrees")

    day_angle = (2.0 * pi * (day_number - 1.0)) / 365.0  # this is in radians
    hour_angle_degrees = 15.0 * (hours_since_midnight + equation_of_time(day_angle) / 60.0 +
                                 ((int(longitude_degrees / 15.0) * 15.0 - longitude_degrees) * 4.0)
                                 / 60.0 + 12.0) - 360.0  # this is in degrees.

    hour_angle = hour_angle_degrees / 180 * numpy.pi  # convert(hour_angle_degrees, "degrees", "radians")

    declination = (
        0.006918
        - 0.399912 * cos(day_angle)
        + 0.070257 * sin(day_angle)
        - 0.006758 * cos(2 * day_angle)
        + 0.000907 * sin(2 * day_angle)
        - 0.002697 * cos(3 * day_angle)
        + 0.00148 * sin(3 * day_angle)
    )

    earth_sun_distance_factor = (
        1.00011
        + 0.034221 * cos(day_angle)
        + 0.001280 * sin(day_angle)
        + 0.000719 * cos(2 * day_angle)
        + 0.000077 * sin(2 * day_angle)
    )

    solar_zenith_angle = arccos(
        cos(declination) * cos(latitude) * cos(hour_angle)
        + sin(declination) * sin(latitude)
    )

    solar_zenith_angle_degrees = solar_zenith_angle / pi * 180  # convert(solar_zenith_angle,"radians",u'degrees')
    # original code checked to stop sun dropping below horizon
    # //Stop sun dropping below horizon
    #     if(solar_zenith_angle > 91)
    #        solar_zenith_angle = 91; (was in degrees)


    # this used to be 1/ could this cause issues in java?
    relative_am = 1.0 / (cos(solar_zenith_angle) + (
        0.15 * pow((93.885 - solar_zenith_angle_degrees), -1.253)))  ##AARG raised to the power of degrees
    pressure_corrected_am = relative_am * (pressure / 101325.33538686013)  # si("1 atm")
    effective_ozone_am = (1.0 + (22.0 / 6370.0)) / (
        pow((pow(cos(solar_zenith_angle), 2.0) + (2.0 * (22.0 / 6370.0))), 0.5))

    if aod_model == "rural":
        c_coefficient = [0.581, 16.823, 17.539]
        d_coefficient = [0.8547, 78.696, 0, 64.458]
    elif aod_model == "rural":
        c_coefficient = [0.2595, 33.843, 39.524]
        d_coefficient = [1.0, 84.254, -9.1, 65.458]
    elif aod_model == "maritime":
        c_coefficient = [0.1134, 0.8941, 1.0796]
        d_coefficient = [0.04435, 1.6048, 0, 1.5298];
    elif aod_model == "tropospheric":
        c_coefficient = [0.6786, 13.899, 13.313]
        d_coefficient = [1.8379, 14.912, 0, 5.96]
    elif type(aod_model) == list:
        c_coefficient, d_coefficient = aod_model

    # 0.9 degrees * something_in_percent rather than 90 degrees?? Honestly.
    x_humidity = cos(humidity * pi / 2.)  # si(0.9,"degrees","radians")

    alpha1 = (c_coefficient[0] + c_coefficient[1] * x_humidity) / (1 + c_coefficient[2] * x_humidity)
    alpha2 = (d_coefficient[0] + d_coefficient[1] * x_humidity + d_coefficient[2] * x_humidity ** 2) / (
        1 + (d_coefficient[3] * x_humidity))

    wavelength_um = am_zero_wavelength * 1e6  # convert(am_zero_wavelength, "m", 'um')

    rayleigh_coeff = exp(-pressure_corrected_am / (wavelength_um ** 4 * (115.6406 - 1.335 / wavelength_um ** 2)))

    aerosol_coeff = where(wavelength_um <= 0.5,
                          exp(-pressure_corrected_am * turbidity * 2.0 ** (alpha2 - alpha1) * wavelength_um ** -alpha1),
                          exp(-pressure_corrected_am * turbidity * wavelength_um ** -alpha2)
                          )

    vapour_coeff = exp(-(0.2385 * precipwater * waterspectra * relative_am)
                       / (1. + 20.07 * precipwater * waterspectra * relative_am) ** 0.45)

    ozone_coeff = exp(-ozonespectra * ozone * effective_ozone_am)

    mixed_coeff = exp((-1.41 * uniformgasspectra * pressure_corrected_am)
                      / (1. + 118.93 * uniformgasspectra * pressure_corrected_am) ** 0.45);

    # Apply attenuation/scatting coefficients to the per m specturm (SI wavelength spectrum)
    wavelength_m = am_zero_wavelength
    irradiance_per_m = am_zero_irradiance * earth_sun_distance_factor * rayleigh_coeff * aerosol_coeff * vapour_coeff * ozone_coeff * mixed_coeff

    if suppress_nan:
        irradiance_per_m = numpy.nan_to_num(irradiance_per_m)

    storage = dict()

    # Convert to per nm spectrum
    wavelength_nm = wavelength_m * 1e9  # asUnit(wavelength_m,"nm")
    irradiance_per_nm = 1e-9 * irradiance_per_m  # asUnit(irradiance_per_m,"nm-1")

    if power_density_in_nm:
        storage = wavelength_nm, irradiance_per_nm

    else:
        # Convert to per eV spectrum
        energy_ev, irradiance_per_ev = spectral_conversion_nm_ev(wavelength_nm, irradiance_per_nm)

        # Covert to energy per Joule spectrum
        energy_j = energy_ev * 1.6e-19  # siUnits(energy_ev, 'J')
        irradiance_per_j = irradiance_per_ev / 1.6e-19  # siUnits(irradiance_per_ev, 'J-1')

        # integrate to find power density
        power_density = trapz(x=wavelength_m, y=irradiance_per_m)

        storage['incident power density'] = power_density
        storage['incident spectrum wavelength si'] = array([wavelength_m, irradiance_per_m])
        storage['incident spectrum wavelength nm'] = array([wavelength_nm, irradiance_per_nm])
        storage['incident spectrum energy si'] = array([energy_j, irradiance_per_j])
        storage['incident spectrum energy eV'] = array([energy_ev, irradiance_per_ev])

    return storage


if __name__ == "__main__":
    data = calculate_spectrum_spectral2()

    import matplotlib.pyplot as plt

    plt.plot(data['incident spectrum wavelength nm'][0], data['incident spectrum wavelength nm'][1])
    plt.show()
