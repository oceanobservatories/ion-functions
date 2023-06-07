import importlib
import requests

# URL to get data from the system.
# since we are not specifying "parameters", all raw data and derived (calculated) data products will be returned.
# if limit == -1 or limit > 1000, request will be serviced asyncronously. Otherwise it is handled syncronously as it is below.
url = "https://ooinet.oceanobservatories.org/api/m2m/12576/sensor/inv/CE02SHBP/LJ01D/05-ADCPTB104/streamed/adcp_velocity_beam?beginDT=2023-06-01T00:00:00.000Z&endDT=2023-06-01T00:01:00.000Z&limit=5"
payload = ""
headers = {'authorization': 'Basic T09JQVBJLVpLV0IyT1RSSjFKR0c4OkYzT1U1V1hMN0tCWjJW'}

# get the data to pass to the ion function
response = requests.request("GET", url, data=payload, headers=headers)

print(response.text)

# convert to python obect
particles = response.json()

# specify the python module (file) and function within the module
pf_owner = "ion_functions.data.adcp_functions"
pf_function = "adcp_bin_depths_meters"

# load the module
module = importlib.import_module(pf_owner)

# stream engine has already called the same function and has returned to us the data product "bin_depths" as well as the raw data inputs
# print it here just to compare it to the results that we will generate just to see that they are the same
print("\nbin_depths as calculated by stream_engine:")
for particle in particles:
    print(particle["bin_depths"])


print("\nbin_depths as calculated here using singular args:")
for particle in particles:
    kwargs = {}
    kwargs["sensor_depth"] = particle["non_zero_depth"]
    kwargs["bin_size"] = particle["cell_length"]
    kwargs["dist_first_bin"] = particle["bin_1_distance"]
    kwargs["num_bins"] = particle["num_cells"]
    kwargs["adcp_orientation"] = particle["sysconfig_vertical_orientation"]

    result = getattr(module, pf_function)(**kwargs)
    print(result)


print("\nbin_depths as calculated here using list args:")
list_kwargs = {}
for particle in particles:
    list_kwargs.setdefault("sensor_depth", []).append(particle["non_zero_depth"])
    list_kwargs.setdefault("bin_size", []).append(particle["cell_length"])
    list_kwargs.setdefault("dist_first_bin", []).append(particle["bin_1_distance"])
    list_kwargs.setdefault("num_bins", []).append(particle["num_cells"])
    list_kwargs.setdefault("adcp_orientation", []).append(particle["sysconfig_vertical_orientation"])

result = getattr(module, pf_function)(**list_kwargs)
print(result)
