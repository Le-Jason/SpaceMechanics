

import requests
import json
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


url = "http://api.open-notify.org/iss-now.json"


response = requests.get(url)


data = json.loads(response.text)


iss_lat = float(data["iss_position"]["latitude"])
iss_lon = float(data["iss_position"]["longitude"])


map = Basemap(projection="cyl")


map.drawcoastlines(linewidth=0.25)
map.fillcontinents(color="coral", lake_color="aqua")


x, y = map(iss_lon, iss_lat)
map.plot(x, y, "bo", markersize=6)


plt.show()