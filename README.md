# PVLIB_PV-System-Installation-for-Arsi-Nagele-General-Hospital

# PV System Installation for Arsi Nagele General Hospital

## Objective
The objective of this project is to install a photovoltaic (PV) system on the rooftop of Arsi Nagele General Hospital to generate solar energy. The south-facing side of the roof was selected for installation.

## Site Assessment
- **Location:** Arsi Nagele General Hospital, Arsi Nagele
- **Methodology:**
  - **Geospatial Analysis:** Google Earth imagery was utilized to extract geospatial data and estimate the rooftop area available for PV module installation. The data was exported as a Keyhole Markup Language (KML) file for further processing.
  - **Meteorological Data Processing:** ERA5 reanalysis climate data was retrieved and processed to obtain historical irradiance, temperature, and weather conditions relevant to the project site.
  - **Computational Data Analysis:** Python scripting was employed to analyze and preprocess the collected spatial and meteorological datasets. The available rooftop area was determined by filtering out obstructions and considering structural constraints.
  - **PV System Simulation:** Once the optimal installation area was identified, PVLiB was utilized to simulate the performance of the PV system under real-world conditions, incorporating parameters such as irradiance, temperature coefficients, and system losses.
- **Roof Orientation:** South-facing

## PV Module Selection & Layout
- **Module Dimensions:** 1.7m x 1m
- **Spacing Between Modules:** 0.1m
- **Total Modules Installed:** 80
- **Layout Determination:** The estimated rooftop area was used to design the module placement, ensuring adequate spacing for ventilation and maintenance. The final configuration was determined through computational analysis, optimizing the number of modules that could be accommodated while maintaining structural integrity.

## System Capacity
- **Nominal DC Power per Module:** 200 W
- **Total System Capacity:** 16 kW (80 modules × 200 W)
- **Azimuth Angle:** 180 degrees (South-facing)

## Simulation Results
The PV system was simulated using PVLiB, and the results were as follows:
- **System Capacity:** 16.00 kW DC
- **Capacity Factor:** 25.05%
- **Specific Yield:** 6.02 kWh/kW
- **Total Energy Generated:** 769.6 kWh (for 8 simulation days)

## Financial Analysis
A basic financial analysis was conducted with the following assumptions:
- **Cost per Watt:** €1.5 EUR/W
- **System Cost:** €24,000.00 EUR

Given that the Levelized Cost of Energy (LCOE) for solar PV is higher than the current electricity tariff in Ethiopia, installing a PV system may not lead to immediate cost savings per kWh for consumers paying standard rates. However, it is crucial to consider the unreliability and frequent interruptions of the Ethiopian electricity grid. For businesses and essential services whose operations depend on a stable and continuous power supply, investing in an off-grid solar PV system can provide energy security, operational stability, and long-term resilience, making it a strategic investment despite the higher initial costs.

## Future Improvements
The model can be enhanced by:
- **Integrating Battery Storage:** To ensure energy availability during non-sunny hours.
- **Demand Profile Optimization:** Aligning system output with the hospital's electricity demand.
- **Cost Analysis:** Comparing grid dependency, fully off-grid systems, and hybrid battery-grid solutions.


