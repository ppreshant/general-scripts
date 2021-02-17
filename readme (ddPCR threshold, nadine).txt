On 31/8/20 7:35 am, Nadine Kotlarz wrote:

Prashant,
I think you would want to use the distribution of negative droplets in a true negative (a sample you know to be negative for the virus). Ideally we would have a matrix matched negative but if you don't have that, maybe you just use the no template control (H2O)? Here are the steps to output the file with fluorescence for each droplet for each sample (must be done on the lab computer itself):

    Open the data in Quantasoft (not AP)
    Go to Setup tab (should be selected by default if the data had to be loaded)
    Select each well individually 
    Select Options (top center)
    Select "Export Amplitude and Cluster Data" and select the folder to save the .csv file


Attached is a script that I asked our biostatistician to write up which calculates threshold using two approaches 1) upper confidence limit around the mean of the negative, 2) 95th percentile
I haven't looked at it closely yet; I'm sharing it in case you find it helpful.

Nadine