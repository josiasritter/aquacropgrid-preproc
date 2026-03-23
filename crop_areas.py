def crop_areas(domain_path, spam_variable, start_year, end_year, basepath, to_match):

    import os
    from preproc_tools import spam_refyear, preproc_spam, makedirs, download_url, unzip_all

    # Define most suitable reference year of crop mask
    refyear = spam_refyear(start_year, end_year)

    ## Download global crop masks from SPAM data (cultivated areas for 16 crop types)
    def download_spam(refyear, spam_variable, basepath):
        # Prepare download URL based on refyear and variable
        if refyear == '2010':
            if spam_variable == 'physical_area':
                url = "https://s3.amazonaws.com/mapspam-data/2010/v2.0/geotiff/spam2010v2r0_global_phys_area.geotiff.zip"
            elif spam_variable == 'yield':
                url = "https://s3.amazonaws.com/mapspam/2010/v2.0/geotiff/spam2010v2r0_global_yield.geotiff.zip"
        elif refyear == '2020':
            if spam_variable == 'physical_area':
                url = "https://www.dropbox.com/scl/fi/napqtql4521ujqt22j05w/spam2020V1r0_global_physical_area.geotiff.zip?rlkey=vpamm4zj3gu2752ubpj3j80iu&e=1&dl=1"
            elif spam_variable == 'yield':
                url = "https://www.dropbox.com/scl/fi/kajp48kh5wnh65ar2ltbr/spam2020V2r0_global_yield.geotiff.zip?rlkey=n1w5823k0ra9uqqg1tbc18ag4&e=1&dl=1"

        # Define target dir and paths
        if spam_variable in ['physical_area', 'harvested_area']:
            target_dir = makedirs(basepath, 'rawdata', 'cropmasks')
        else:  # yield, production
            target_dir = makedirs(basepath, 'rawdata', 'calibration')

        download_path = os.path.join(target_dir, f'spam{refyear}_{spam_variable}.zip')
        unzipped_dir = download_path[:-4]

        # If unzipped data already exists, skip everything
        if os.path.exists(unzipped_dir):
            print(f"SPAM data already unzipped, skipping download and unzip: {unzipped_dir}")
            return unzipped_dir

        # Otherwise, check if ZIP exists
        if os.path.exists(download_path):
            print(f" SPAM zip already exists, skipping download: {download_path}")
        else:
            print(f"  Downloading SPAM {refyear} data ({spam_variable})")
            print('    URL:', url)
            download_url(url, download_path=download_path)

        # Unzip (if not already unzipped)
        print(" Unzipping SPAM data...")
        unzip_all(dir=target_dir)

        return unzipped_dir

    # Run downloader
    download_dir = download_spam(refyear, spam_variable, basepath)

    ## Preprocess data for model domain
    preproc_spam(basepath, download_dir, refyear, spam_variable, domain_path, to_match)
