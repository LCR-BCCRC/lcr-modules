import sys
import pkg_resources


def checker():
# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
    min_oncopipe_version="1.0.11"
    print("ss")
    try:
        from packaging import version
    except ModuleNotFoundError:
        sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

    current_version = pkg_resources.get_distribution("oncopipe").version
    if version.parse(current_version) < version.parse(min_oncopipe_version):
        print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
        print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
        sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section 

checker()