import os
import pooch
import zipfile
import shutil


def get_zenodo(url, save_path, checksum, processor=None):

    print(f"\nChecking if {save_path} is already downloaded...")

    # Ensure the save directory exists
    os.makedirs(os.path.dirname(save_path), exist_ok=True)

    # Fetch the file with pooch
    file_path = pooch.retrieve(
        url=url,
        path=os.path.split(save_path)[0],
        # Using a cache directory under pooch's management
        known_hash=f"md5:{checksum}",
        fname=os.path.basename(save_path),
        progressbar=True,
        processor=processor,
    )


def download_data():
    fair_calibrate = "datasets/fair_calibrate"

    get_zenodo(
        "https://zenodo.org/records/10566813/files/fair_calibrate.zip",
        os.path.join(fair_calibrate, "fair_calibrate.zip"),
        "2509a7a03b81f84a19d4e3aa7f9e65fe",
        pooch.Unzip(extract_dir="."),
    )

    # RCMIP = "dependencies/RCMIP"

    # get_zenodo(
    #     "https://zenodo.org/records/4589756/files/rcmip-emissions-annual-means-v5-1-0.csv",
    #     os.path.join(RCMIP, "rcmip-emissions-annual-means-v5-1-0.csv"),
    #     "4044106f55ca65b094670e7577eaf9b3",
    # )

    # get_zenodo(
    #     "https://zenodo.org/records/4589756/files/rcmip-concentrations-annual-means-v5-1-0.csv",                                                                          
    #     os.path.join(RCMIP, "rcmip-concentrations-annual-means-v5-1-0.csv"), 
    #     "0d82c3c3cdd4dd632b2bb9449a5c315f",
    # )

    # get_zenodo(
    #     "https://zenodo.org/records/4589756/files/rcmip-radiative-forcing-annual-means-v5-1-0.csv",
    #     os.path.join(RCMIP, "rcmip-radiative-forcing-annual-means-v5-1-0.csv"), 
    #     "87ef6cd4e12ae0b331f516ea7f82ccba",
    # )


if __name__ == "__main__":
    download_data()
