import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, quantity_support
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon, name_resolve


def pltBody(body, midnight, loc, delta_settime):
    # Find altaz in set timeframe
    frame_night = AltAz(
        obstime=midnight + delta_settime, location=loc)

    # Set body coordinates
    objectaltazs_night = body.transform_to(frame_night)

    # Convert alt, az to airmass using secz
    objectairmass_night = objectaltazs_night.secz  # secz estimates airmass

    # Plot airmass as a function of time, using matplotlib
    plt.plot(delta_settime, objectairmass_night)
    plt.xlim(9, 23)
    plt.ylim(1, 8)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Airmass [Sec(z)]')
    plt.show()


def pltSMB(body, bodyName, midnight, loc, delta_settime):
    '''
        Plots the relationship between the Sun / Moon and the specified body
    '''

    # Find the Sun at 1000 spaced times between noon July 12 and noon July 13
    delta_midnight = np.linspace(-12, 12, 1000) * u.hour
    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(
        obstime=times_July12_to_13, location=loc)
    sunaltazs_July12_to_13 = get_sun(
        times_July12_to_13).transform_to(frame_July12_to_13)

    # Do the same with get_moon
    moon_July12_to_13 = get_moon(times_July12_to_13)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(
        frame_July12_to_13)

    # Find alt,az coordinates of body at these times
    tonaltazs_July12_to_13 = body.transform_to(frame_July12_to_13)

    # Use matplotlib to make a figure illustrating nighttime and altitudes
    plt.plot(delta_midnight, sunaltazs_July12_to_13.alt, color='r', label='Sun')
    plt.plot(delta_midnight, moonaltazs_July12_to_13.alt,
             color=[0.75] * 3, ls='--', label='Moon')
    plt.scatter(delta_midnight, tonaltazs_July12_to_13.alt,
                c=tonaltazs_July12_to_13.az, label=bodyName, lw=0, s=8, cmap='viridis')

    plt.fill_between(delta_midnight, 0 * u.deg, 90 * u.deg,
                     sunaltazs_July12_to_13.alt < -18 * u.deg, color='k', zorder=0)
    plt.colorbar().set_label('Azimuth [deg]')
    plt.legend(loc='upper left')
    plt.xlim(-12 * u.hour, 12 * u.hour)
    plt.xticks((np.arange(13) * 2 - 12) * u.hour)
    plt.ylim(0 * u.deg, 90 * u.deg)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Altitude [deg]')
    plt.show()


def main():
    plt.style.use(astropy_mpl_style)
    quantity_support()

    bodyName = " "
    choice = 0

    while choice != 3 and bodyName != "":
        bodyName = input("Enter celestial body name (nothing to exit): ")

        try:
            body = SkyCoord.from_name(bodyName)
        except name_resolve.NameResolveError:
            if bodyName != "":
                print("E: Celestial body not found.")
        else:
            if body != "":
                # Get location coordinates from the user
                print("Please enter the observation point location (latitude/longitude/height: )")
                latitude = input("Lat: ")
                longitude = input("Long: ")
                height = input("Height(m): ")

                # Convert location coordinates to floats
                try:
                    latitude = float(latitude)
                    longitude = float(longitude)
                    height = int(height)
                except ValueError:
                    print("Latitude/Longitude/Height values not correct.")
                else:
                    # Set location of observation using astropy.units
                    location = EarthLocation(
                        lat=latitude * u.deg, lon=-longitude * u.deg, height=height * u.m)

                    # lat=41.3 * u.deg, lon=-74 * u.deg, height=390 * u.m

                    utcoffset = -4 * u.hour  # Eastern Daylight Time

                    # User input for date and time
                    date = input("Enter date (YYY-M-D): ")
                    time = input("Enter time in EDT (HH:MM:SS): ")

                    datetime = Time(date + " " + time) - \
                        utcoffset  # Set observation time

                    # Find coordinates of selected body as observed at set time
                    tonaltaz = body.transform_to(
                        AltAz(obstime=datetime, location=location))
                    print(f"{bodyName}'s Altitude = {tonaltaz.alt:.6}")

                    midnight = Time(date + ' 00:00:00') - utcoffset
                    # linspace returns evenly spaced numbers over specified interval
                    delta_settime = np.linspace(9, 23, 100) * u.hour

                    choice = input(
                        "Plot Body (1) | Plot Moon & Sun vs Body (2) | Quit (3): ")
                    try:
                        choice = int(choice)
                    except ValueError:
                        print("Please enter a valid option.")
                    else:
                        if choice == 1:
                            pltBody(body, midnight, location, delta_settime)
                        elif choice == 2:
                            pltSMB(body, bodyName, midnight,
                                   location, delta_settime)
                        else:
                            if choice != 3:
                                print("Please enter a valid option.")


# Call main function
if __name__ == "__main__":
    main()
