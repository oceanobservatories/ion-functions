import importlib
import requests

# URL to get data from the system.
# since we are not specifying "parameters", all raw data and derived (calculated) data products will be returned.
# if limit == -1 or limit > 1000, request will be serviced asyncronously. Otherwise it is handled syncronously as it is below.
url = "https://ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv/CE01ISSM/RID16/03-CTDBPC000/telemetered/ctdbp_cdef_dcl_instrument?beginDT=2023-06-01T00:00:00.000Z&endDT=2023-06-02T00:00:00.000Z&limit=5"
payload = ""
headers = {'authorization': 'Basic T09JQVBJLVpLV0IyT1RSSjFKR0c4OkYzT1U1V1hMN0tCWjJW'}

# get the data to pass to the ion function
response = requests.request("GET", url, data=payload, headers=headers)

print(response.text)

# convert to python obect
particles = response.json()

# specify the python module (file) and function within the module
pf_owner = "ion_functions.data.ctd_functions"
pf_function = "ctd_pracsal"

# load the module
module = importlib.import_module(pf_owner)

# stream engine has already called the same function and has returned to us the data product "bin_depths" as well as the raw data inputs
# print it here just to compare it to the results that we will generate just to see that they are the same
print("\nsalinity as calculated by stream_engine:")
for particle in particles:
    print(particle["sea_water_practical_salinity"])


        # c = sea water conductivity (CONDWAT_L1) [S m-1]
        # t = sea water temperature (TEMPWAT_L1) [deg_C]
        # p = sea water pressure (PRESWAT_L1) [dbar]
print("\nsalinity as calculated here using singular args:")
for particle in particles:
    kwargs = {}
    kwargs["c"] = particle["sea_water_electrical_conductivity"]
    kwargs["t"] = particle["sea_water_temperature"]
    kwargs["p"] = particle["sea_water_pressure"]

    result = getattr(module, pf_function)(**kwargs)
    print(result)


# print("\nsalinity as calculated here using list args:")
# list_kwargs = {}
# for particle in particles:
#     list_kwargs.setdefault("c", []).append(particle["sea_water_electrical_conductivity"])
#     list_kwargs.setdefault("t", []).append(particle["sea_water_temperature"])
#     list_kwargs.setdefault("p", []).append(particle["sea_water_pressure"])
#
# result = getattr(module, pf_function)(**list_kwargs)
# print(result)
