def_feat_tab_fname <- "Features_aligned_merged_IceR_analysis.txt"

#' Function to determine at which value a density maximum is reached
#' @param data Numeric vector
#' @details Uses kernel density estimation function from R-package stats
#' removing missing values
#' @return Numeric indicating at which value a density maximum is reached
#' @export
maxDensity <- function(data){
    dens <- stats::density(data, na.rm=TRUE)
    return(dens$x[which(dens$y == max(dens$y))])
}

#' Convert thermo raw files to mzXML with centroided ms1 scans using the
#' ProteoWizard tool msConvert
#' @param path_to_raw Path to folder containing raw files which should be
#' converted
#' @details Requires installation of ProteoWizard
#' '(http://proteowizard.sourceforge.net/download.html). Pay attention to
#' installation requirements.
#' @return Resulting mzXML files are stored in a sub-directory within specified
#' raw file folder
#' @export
run_msconvert_raw_mzXML <- function(path_to_raw=NULL){
    if(is.null(path_to_raw)){
        path_to_raw <- rChoiceDialogs::rchoose.dir(
            caption="Select folder containing raw files"
        )
    }

    raw_files <- list.files(path_to_raw)
    raw_files <- raw_files[which(grepl("\\.raw", raw_files))]

    print("The following raw files were found in the specified folder: ")
    cat(raw_files, sep='\n')

    # Check which raw files still have to be converted
    mzXMLs_available <- list.files(base::paste(path_to_raw, "/mzXML", sep=""))

    c1 <- which(base::gsub("\\.raw", "", raw_files)
    c2 <- base::gsub("\\.mzXML", "", mzXMLs_available))
    files_to_be_converted <- raw_files[c1 %not in% c2]

    if(length(files_to_be_converted) > 0){
        # Get home directory
        home_folder <- Sys.getenv("HOME")
        home_folder <- base::gsub("Documents", "AppData/Local/Apps/",
                                  home_folder)

        # Find MSConvert folder
        folders <- list.dirs(path=home_folder, full.names=TRUE, recursive=FALSE)
        folders <- folders[which(grepl("ProteoWizard", folders))]
        folders <- folders[length(folders)]

        if(file.exists(base::paste(folders, "\\msconvert.exe", sep=""))){
            path_to_msconvert <- base::paste(folders, "\\msconvert.exe", sep="")
        }

        else{
            pth <- "C:/Users/user_name/AppData/Local/Apps/ProteoWizard Version"
            print(paste0("Could not find msconvert.exe. Can be usually found ",
                         "in ", pth))
            pb <- tcltk::tkProgressBar("Warning!", min=0, max=3, initial=0,
                                       label="Could not find msConvert.exe.",
                                       width=500)

            counter <- 1

            label <- c(paste0("On Windows typically located in ", pth),
                       "Please specify location of msConvert.exe")

            while(TRUE){
                Sys.sleep(3)
                tcltk::setTkProgressBar(pb, value=counter, label=label[counter])
                counter <- counter + 1

                if(counter == 4){
                    break
                }
            }

            close(pb)
            path_to_msconvert <- file.choose()
        }

        # Get user folder
        win_user_folder <- path.expand('~')

        # Create temporary folder
        dir.create(base::paste(win_user_folder, "\\temp_msconvert", sep=""),
                   showWarnings=FALSE)
        dir.create(base::paste(win_user_folder, "\\temp_msconvert\\mzXML",
                               sep=""), showWarnings=FALSE)

        # Create temporary config.txt and files.txt
        temp_path <- base::paste(win_user_folder, "\\temp_msconvert", sep="")
        setwd(temp_path)
        # Config
        fileConn <- file("config.txt")
        writeLines(c("mzXML=true", "64=true", "noindex=false", "zlib=true",
                     "filter=\"peakPicking vendor msLevel=1\"",
                     "filter=\"msLevel 1\""), fileConn)
        close(fileConn)
        # Files
        fileConn <- file("files.txt")
        writeLines(base::paste(path_to_raw, "\\", files_to_be_converted,
                               sep=""), fileConn)
        close(fileConn)

        # Prepare arguments for msconvert
        arg <- base::paste("-f ", temp_path, "\\files.txt", " -o ", temp_path,
                           "\\mzXML", " -c ", temp_path, "\\config.txt", sep="")
        # Run msconvert
        system2(path_to_msconvert, args=arg)

        # Move mzXMLs to original raw folder
        from <- temp_path
        to <- path_to_raw
        path1 <- base::paste0(from, "\\mzXML")
        path2 <- base::paste0(to, "\\mzXML")

        dir.create(path2, showWarnings=FALSE)

        for(f in base::gsub("\\.raw", ".mzXML", files_to_be_converted)){
            ff::file.move(base::paste(path1, "\\", f, sep=""), path2)
        }

        # Remove temp msconvert folder
        r <- file.remove(c("files.txt", "config.txt"))
    }
}

#' Prepare mzXML files for IceR workflow
#' @param path_to_mzXML Path to folder containing mzXML files
#' @param n_cores Numbers of CPU cores which should be used to perform
#' conversion
#' @details Converts MS1 spectra in mzXML files into tables containing m/z, RT
#' and intensity information per ion
#' @return Resulting ion tables are stored in a sub-directory (all_ion_lists) of
#' the mzXML folder as .RData files
#' @export
mzxml_to_list <- function(path_to_mzXML, n_cores=2){
    i <- 0

    convert <- function(mzXMLfile, path_to_mzXML){
        data <- base::paste(path_to_mzXML, "/", mzXMLfile, sep="")
        label <- base::paste(round(0 / 1 * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title="Read mzXML", label=label, min=0,
                                   max=1, width=300)
        ms <- readMzXmlData::readMzXmlFile(data)
        close(pb)

        # Get total rowcount
        rowcount <- 0

        for(i in 1:length(ms)){
            if(ms[[i]]$metaData$msLevel == 1){
                rowcount <- rowcount + length(ms[[i]]$spectrum$mass)
            }
        }

        # Now extract data
        dat <- data.table::as.data.table(matrix(ncol=3, nrow=rowcount))
        colnames(dat) <- c("m.z", "RT", "Intensity")
        sample <- mzXMLfile
        sample <- base::substr(sample, 1, regexpr(".mzXML", sample) - 1)
        dat$m.z <- as.numeric(dat$m.z)
        dat$RT <- as.numeric(dat$RT)
        dat$Intensity <- as.numeric(dat$Intensity)
        ind <- 1
        max <- length(ms)
        label <- base::paste(round(0 / max * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title="Extract data", label=label, min=0,
                                   max=max, width=300)

        for(i in 1:length(ms)){
            if(ms[[i]]$metaData$msLevel == 1){
                start <- ind
                stop <- ind + length(ms[[i]]$spectrum$mass) - 1
                data.table::set(x=dat, i=start:stop, j=(1L),
                                value=ms[[i]]$spectrum$mass)
                data.table::set(x=dat, i=start:stop, j=(2L),
                                value=ms[[i]]$metaData$retentionTime / 60)
                data.table::set(x=dat, i=start:stop, j=(3L),
                                value=ms[[i]]$spectrum$intensity)
                ind <- ind + length(ms[[i]]$spectrum$mass)
            }

            label <- base::paste(round(i / max * 100, 0), "% done (", i, "/",
                                 max, ")", sep="")
            tcltk::setTkProgressBar(pb, i, label=label)
        }

        close(pb)
        save(dat, file=base::paste(path_to_mzXML, "/all_ion_lists/", sample,
                                   "_all_ions.RData", sep=""))
    }

    ##Step - Extract all ions per ms1 spectra

    mzXMLfiles <- list.files(path_to_mzXML)
    mzXMLfiles <- mzXMLfiles[which(grepl(".mzXML", mzXMLfiles))]

    dir.create(base::paste(path_to_mzXML, "/all_ion_lists", sep=""),
               showWarnings=FALSE)

    setwd(base::paste(path_to_mzXML, "/all_ion_lists", sep=""))

    # Check if all ion list are already available and only generate those which
    # are required
    loc <- which(grepl("\\.RData", list.files()))
    available_all_ion_lists <- list.files()[loc]

    c1 <- base::gsub("\\.mzXML", "", mzXMLfiles)
    c2 <- base::gsub("_all_ions\\.RData", "", available_all_ion_lists)
    missing_all_ion_lists <- mzXMLfiles[which(c1 %not in% c2)]
    if(length(missing_all_ion_lists) > 0){
        mzXMLfiles <- missing_all_ion_lists

        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)
        res <- foreach::foreach(i=mzXMLfiles) %dopar% {
            convert(i, path_to_mzXML)
        }

        parallel::stopCluster(cl)
    }
}

#' Extract MS1-spectra from raw TIMS-ToF Pro data
#' @param path_to_raw Path to folder containing d-files
#' @details Convert MS1 spectra with TIMS information in raw files of Bruker
#' TIMS-ToF Pro Mass Spectrometers into tables containing m/z, RT, inverse ion
#' mobility and intensity information per ion. This process is currently very
#' slow and can take several days even for a small data set. It also requires
#' enough space on the home drive of the PC as well as at least 50 Gb of memory.
#' @return Resulting ion tables are stored in a sub-directory (all_ion_lists) of
#' the raw folder as .RData files
#' @export
convert_rawTIMS <- function(path_to_raw=NULL){
    # Extract spectra from raw text
    extract_spectra <- function(path, filename){
        con <- file(base::paste(path, filename, sep=""), "r")

        # Read until spectrumList length is given
        while(TRUE){
            line <- readLines(con, 1)

            if(grepl("    spectrumList \\(", line)){
                break
            }
        }

        # Get number of total spectra
        n_spectra <- as.numeric(base::gsub("    spectrumList \\(| spectra):",
                                           "", line[1]))

        # Prepare table summarizing relevant information
        spectra <- base::as.data.frame(matrix(ncol=5, nrow=n_spectra))
        colnames(spectra) <- c("RT", "1/KO", "mz", "int", "num_ions")
        spectra$RT <- as.numeric(spectra$RT)
        spectra$`1/KO` <- as.numeric(spectra$`1/KO`)
        spectra$mz <- as.character(spectra$mz)
        spectra$int <- as.character(spectra$int)
        spectra$num_ions <- as.numeric(spectra$num_ions)
        count_spectra <- 0

        # Jump to section where individual spectra are starting
        while(TRUE){
            line <- readLines(con, 1)

            if(grepl("      spectrum:", line)){
                break
            }
        }

        # Now read individual spectra
        print(base::paste(base::gsub(".txt", "", filename),
                          ": Prepare TIMS-ToF data", sep=""))
        pb <- utils::txtProgressBar(min=0, max=n_spectra, style=3)

        for(i in 1:n_spectra){
            # Check if ion data is available
            dal <- "        defaultArrayLength: "

            if(i == 1){
                line <- readLines(con, 3)
                ion_count <- as.numeric(base::gsub(dal, "", line[3]))
            }

            else{
                line <- readLines(con, 4)
                ion_count <- as.numeric(base::gsub(dal, "", line[4]))
            }

            if(ion_count > 0){
                # Read data for this spectrum
                line <- readLines(con, 43)

                # Increment number of available spectra with ions
                count_spectra <- count_spectra + 1

                # Get RT and 1/K0
                x <- "            cvParam: scan start time, |, second"
                RT <- as.numeric(base::gsub(x, "", line[32]))
                x <- paste0("            cvParam: inverse reduced ion ",
                            "mobility, |, volt-second per square centimeter")
                ims <- as.numeric(base::gsub(x, "", line[33]))

                # Get m/z and intensities
                mz <- base::gsub(base::paste("          binary: \\[",
                                 ion_count, "\\] ", sep=""), "", line[40])
                int <- base::gsub(base::paste("          binary: \\[",
                                  ion_count, "\\] ", sep=""), "", line[43])

                # Store data in table
                data.table::set(spectra, as.integer(count_spectra),
                                as.integer(1:5), list(RT, ims, mz, int,
                                                      ion_count))
            }

            else{
                # Jump to the end of this spectrum
                line <- readLines(con, 41)
            }

            utils::setTxtProgressBar(pb, i)
        }

        close(pb)
        spectra <- spectra[1:count_spectra, ]
        close(con)

        # Prepare final table
        table_store <- base::as.data.frame(matrix(ncol=4,
                                                  nrow=sum(spectra$num_ions)))
        colnames(table_store) <- c("RT", "1/KO", "mz", "int")
        table_store$RT <- as.numeric(table_store$RT)
        table_store$`1/KO` <- as.numeric(table_store$`1/KO`)
        table_store$mz <- as.numeric(table_store$mz)
        table_store$int <- as.numeric(table_store$int)
        counter <- 1
        print(base::paste(base::gsub(".txt", "", filename),
                          ": Finalize conversion of TIMS-ToF data", sep=""))
        pb <- utils::txtProgressBar(min=0, max=nrow(spectra), style=3)

        for(i in 1:nrow(spectra)){
            start <- counter
            end <- counter + spectra$num_ions[i] - 1
            mzs <- as.numeric(unlist(strsplit(spectra$mz[i], " ")))
            ints <- as.numeric(unlist(strsplit(spectra$int[i], " ")))
            data.table::set(table_store, as.integer(start:end), as.integer(1:4),
                            list(rep(spectra$RT[i], spectra$num_ions[i]),
                                 rep(spectra$`1/KO`[i], spectra$num_ions[i]),
                                 mzs, ints))
            counter <- counter + spectra$num_ions[i]
            utils::setTxtProgressBar(pb, i)
        }

        close(pb)
        dir.create(base::paste(path, "all_ion_lists", sep=""),
                   showWarnings=FALSE)
        print(base::paste(base::gsub(".txt", "", filename),
                          ": Store extracted spectra data", sep=""))
        save(table_store, file=base::paste(path, "all_ion_lists/",
                                           base::gsub(".txt", "_all_ions.RData",
                                                      filename), sep=""))
        rm(spectra,table_store)
        gc()
    }

    '%!in%' <- function(x, y){
        !('%in%'(x, y))
    }

    '%not in%' <- function(x, y){
        !('%in%'(x, y))
    }

    if(is.null(path_to_raw)){
        path_to_raw <- rChoiceDialogs::rchoose.dir(
            caption="Select folder containing raw files"
        )
    }

    # Get home directory
    home_folder <- Sys.getenv("HOME")
    home_folder <- base::gsub("Documents", "AppData/Local/Apps/", home_folder)

    # Find MSConvert folder
    folders <- list.dirs(path=home_folder, full.names=TRUE, recursive=FALSE)
    folders <- folders[which(grepl("ProteoWizard", folders))]
    folders <- folders[length(folders)]

    raw_files <- list.files(path_to_raw)
    raw_files <- raw_files[which(grepl("\\.d", raw_files))]

    # Check which raw files still have to be converted
    text_available <- list.files(base::paste(path_to_raw, "/all_ion_lists",
                                             sep=""))
    c1 <- base::gsub("\\.d", "", raw_files)
    c2 <- base::gsub("_all_ions.RData", "", text_available)
    files_to_be_converted <- raw_files[which(c1 %not in% c2)]

    if(length(files_to_be_converted) > 0){
        # Find location of msconvert.exe
        pth <- "C:/Users/user_name/AppData/Local/Apps/ProteoWizard Version"

        if(file.exists(base::paste(folders, "\\msconvert.exe", sep=""))){
            path_to_msconvert <- base::paste(folders, "\\msconvert.exe", sep="")
        }

        else{
            print(paste0("Could not find msconvert.exe. Can be usually ",
                         "found in ", pth))
            pb <- tcltk::tkProgressBar("Warning!", min=0, max=3, initial=0,
                                       label="Could not find msConvert.exe.",
                                       width=500)

            counter <- 1
            label <- c(paste0("Typically located in ", pth),
                       "Please specify location of msConvert.exe")
            while(TRUE){
                Sys.sleep(3)
                tcltk::setTkProgressBar(pb, value=counter, label=label[counter])
                counter <- counter + 1

                if(counter == 4){
                    break
                }
            }

            close(pb)
            path_to_msconvert <- file.choose()
        }

        # Prepare conversion for each file into raw text
        # Get user folder
        win_user_folder <- path.expand('~')

        # Create temporary folder
        dir.create(base::paste(win_user_folder, "\\temp_msconvert", sep=""),
                   showWarnings=FALSE)
        dir.create(base::paste(win_user_folder, "\\temp_msconvert\\temp",
                               sep=""), showWarnings=FALSE)

        # Create temporary config.txt and files.txt
        temp_path <- base::paste(win_user_folder, "\\temp_msconvert", sep="")
        setwd(temp_path)

        # Config
        fileConn <- file("config.txt")
        writeLines(c("text=true", "64=true", "zlib=true",
                     "filter=\"peakPicking vendor msLevel=1\"",
                     "filter=\"msLevel 1\""), fileConn)
        close(fileConn)

        # Now perform conversion using msconvert and subsequently by R for each
        # file step by step
        for(i in 1:length(files_to_be_converted)){
            # Specify files.txt for msconvert containing the path to the current
            # file
            fileConn <- file("files.txt")
            writeLines(base::paste(path_to_raw, "\\", files_to_be_converted[i],
                                   sep=""), fileConn)
            close(fileConn)

            # Prepare arguments for msconvert
            arg <- base::paste("-f ", temp_path, "\\files.txt", " -o ",
                               temp_path, "\\temp", " -c ", temp_path,
                               "\\config.txt", sep="")

            # Run msconvert
            # Check if file is already converted to raw text
            usedir <- base::paste(temp_path, "\\temp\\",
                                  base::gsub(".d", ".txt",
                                             files_to_be_converted[i]), sep="")

            if(file.exists(usedir) == FALSE){
            # If not use msconvert
                system2(path_to_msconvert, args=arg)
            }

            else{
                print("Converted txt file already available")
            }

            # When finished, trigger conversion into spectra RData file
            # Check if file is already converted to raw text
            usedir <- base::paste(temp_path, "\\temp\\all_ion_lists\\",
                                  base::gsub(".d", "_all_ions.RData",
                                             files_to_be_converted[i]), sep="")

            if(file.exists(usedir) == FALSE){
                # If not use msconvert
                extract_spectra(path=base::paste(temp_path, "\\temp\\", sep=""),
                                filename=base::gsub(".d", ".txt",
                                                    files_to_be_converted[i]))
            }

            else{
                print("Spectra were already extracted")
            }

            # Copy final spectra .RData file to destination and remove raw txt
            # file
            from <- temp_path
            to   <- path_to_raw
            path1 <- base::paste0(from, "\\temp\\all_ion_lists")
            path2 <- base::paste0(to, "\\all_ion_lists")

            dir.create(path2,showWarnings=FALSE)

            f <- base::gsub("\\.d", "_all_ions.RData", files_to_be_converted[i])
            r <- ff::file.move(base::paste(path1, "\\", f, sep=""), path2)
            # Remove raw txt file
            r <- file.remove(base::paste(temp_path, "\\temp\\",
                                         base::gsub(".d", ".txt",
                                                    files_to_be_converted[i]),
                                         sep=""))
        }
    }

    else{
        print("All files are already converted.")
    }
}

#' Perform alignment of pre-determined MS1-features by MaxQuant over proteomics
#' samples either analyzed on an Orbitrap machine (e.g. Q-Exactive or Orbitrap
#' Fusion) or a TIMS-ToF Pro
#' @param path_to_MaxQ_output Path to folder containing MaxQuant outputs (txt
#' folder containing at least allpeptides.txt, evidence.txt, peptides.txt and
#' proteinGroups.txt)
#' @param path_to_output Path to folder where IceR results should be stored
#' @param output_file_names_add IceR result name tag. By default IceR_analysis
#' @param mz_window Numeric value indicating maximal m/z deviation around a
#' center. By default set to NA which indicates that the function determines
#' automatically this parameter based on sd of m/z of identified peptides
#' between samples.
#' @param min_mz_window Numberic value indicating how large the automatically
#' determined m/z-window should be. By default set to 0.001 Da. Only required if
#' m/z alignment window should be automatically determined. If set to NA, no
#' minimal window size will be required.
#' @param RT_window Numeric value indicating maximal RT deviation around a
#' center. By default set to NA which indicates that the function determines
#' automatically this parameter based on sd of RT of identified peptides between
#' samples.
#' @param min_RT_window Numberic value indicating how large the automatically
#' determined RT-window should be. By default set to 1 min. Only required if RT
#' alignment window should be automatically determined. If set to NA, no minimal
#' window size will be required.
#' @param feature_mass_deviation_collapse Numeric value indicating which minimal
#' mass deviation is required to distinguish IceR features. By default set to
#' 0.002 Da. IceR features with overlapping RT-windows and mass differences
#' smaller than the specified value are merged.
#' @param only_unmodified_peptides Boolean value indicating if only unmodified
#' peptide sequences are used for alignment. By default set to F.
#' @param remove_contaminants Boolean value indicating if peptide features
#' labeled as contaminants should be removed. By default set to T.
#' @param sample_list Character vector (raw file names) listing which samples
#' should be aligned. By default all samples occuring in MaxQuant outputs are
#' aligned.
#' @param align_unknown Boolean value indicating if only peptide features or
#' also unsequenced features should be aligend over samples. By default set to
#' FALSE.
#' @param min_num_ions_collapse Numeric value indicating how many unsequenced
#' MaxQuant features have to be at least detected over all samples to result in
#' an IceR feature. Only required if align_unknown is set to T. By default set
#' to 10.
#' @param MassSpec_mode String being either "Orbitrap" or "TIMSToF" specifying
#' by which type of Mass Spectrometer the data was generated. By default it
#' expects Thermo Orbitrap data.
#' @param IM_window Numeric value indicating maximal ion mobility (inverse K0)
#' deviation around a center. By default set to NA which indicates that the
#' function determines automatically this parameter based on sd of ion mobility
#' of identified peptides between samples.
#' @param min_IM_window Numberic value indicating how large the automatically
#' determined ion mobility-window (inverse K0) should be. By default set to
#' 0.002. Only required if ion mobility alignment window should be automatically
#' determined. If set to NA, no minimal window size will be required.
#' @details Performs the first steps of the IceR workflow: 1) Alignment window
#' determination if not specified. 2) Alignment of MaxQuant features into IceR
#' features. 3) Transfer of sequence information between MaxQuant features
#' aligned into IceR features. 4) Extraction, modelling and prediction of RT-
#' and m/z-correction factor per IceR feature and sample.
#' @return Outputs are stored in the sub-directory Temporary_files within
#' specified output folder. MaxQuant allpeptides.txt and evidence.txt are
#' converted to RData files. QC plots of estimated alignment windows as well as
#' of random forest modesl and generalized additive models are stored in a
#' QC_plots.pdf. Relevant QC data is stored in Feature_alignment_QC_data.RData.
#' Aligned IceR features are stored in Features_aligned_merged.txt
#' @export
align_features <- function(path_to_MaxQ_output, path_to_output,
                           align_unknown=FALSE,
                           output_file_names_add="IceR_analysis", mz_window=NA,
                           min_mz_window=0.001, RT_window=NA, min_RT_window=1,
                           min_num_ions_collapse=10,
                           feature_mass_deviation_collapse=0.002,
                           only_unmodified_peptides=FALSE, sample_list=NA,
                           remove_contaminants=TRUE,
                           MassSpec_mode=c("Orbitrap", "TIMSToF"), IM_window=NA,
                           min_IM_window=0.002){
    options(warn=-1)
    use_mz_at_max_int_for_correction <- FALSE

    # Mass spec Mode can only be Orbitrap or TIMSToF
    MassSpec_mode <- MassSpec_mode[1]

    if(is.na(sample_list)){
        sample_list <- NULL
    }

    if(output_file_names_add != ""){
        output_file_names_add <- base::paste("_", output_file_names_add, sep="")
    }

    setwd(path_to_output)
    dir.create("Temporary_files")

    pth <- base::paste("Temporary_files/Features_aligned_merged",
                       output_file_names_add, ".txt", sep="")

    if(file.exists(pth)){
        options(warn=1)
        warning(base::paste("Features_aligned_merged", output_file_names_add,
                            ".txt is already available in the Temporary_files ",
                            "folder\n IceR will use this one for subsequent ",
                            "steps. If this is not intended, please remove ",
                            "Features_aligned_merged", output_file_names_add,
                            ".txt in Temporary_files folder or the complete ",
                            "folder.", sep=""))
        options(warn=-1)
    }

    else{
        if(file.exists("Temporary_files/allPeptides.RData")){
            load("Temporary_files/allPeptides.RData")
            load("Temporary_files/evidence.RData")

            if(is.null(sample_list)){
                sample_list <- sort(unique(allpeptides$Raw.file))
            }

            else{
                # Check if all Raw.files are in sample_list and that all samples
                # in sample_list are in Raw.file
                all_samples_in_MaxQ <- sort(unique(allpeptides$Raw.file))

                # Any samples specified for IceR but not available in MaxQ data?
                if(any(sample_list %not in% all_samples_in_MaxQ)){
                    missing <- which(sample_list %not in% all_samples_in_MaxQ)
                    options(warn=1)
                    miss <- paste(sample_list[missing], collapse=",")
                    warning(paste0(miss, " is/are specified to be ",
                                   "requantified by IceR but missing in ",
                                   "MaxQuant results. This/these samples will ",
                                   "be skipt by IceR."))
                    options(warn=-1)
                }

                # Any samples available in MaxQ but not specified for IceR?
                if(any(all_samples_in_MaxQ %not in% sample_list)){
                    missing <- which(all_samples_in_MaxQ %not in% sample_list)
                    options(warn=1)
                    miss <- paste(all_samples_in_MaxQ[missing], collapse=",")
                    warning(paste0(miss," is/are available in MaxQuant ",
                                   "results but not specified to be ",
                                   "requantified by IceR. This/these samples ",
                                   "will be skipt by IceR."))
                    options(warn=-1)
                }
            }
        }

        else{
            setwd(path_to_MaxQ_output)
            options(fftempdir=path_to_MaxQ_output)
            print(paste0(Sys.time(), " Read MaxQ results"))
            temp <- utils::read.csv(file="allPeptides.txt", sep='\t',
                                    nrows=2044, header=TRUE)
            tempclasses <- sapply(temp, class)

            if(MassSpec_mode == "Orbitrap"){
                tempclasses["Intensity"] <- "numeric"
                tempclasses[which(tempclasses == "logical"
                                  | tempclasses == "character")] <- "factor"
            }

            else{
                tempclasses["m.z"] <- "numeric"
                tempclasses["Mass"] <- "numeric"
                tempclasses["Retention.time" ] <- "numeric"
                tempclasses["Retention.length"] <- "numeric"
                tempclasses["Ion.mobility.index"] <- "numeric"
                tempclasses["Ion.mobility.index.length"] <- "numeric"
                tempclasses[which(tempclasses == "logical"
                                  | tempclasses == "character")] <- "factor"
            }

            # Read in data in chunks of 100000 rows
            allpeptides_save <- ff::read.csv.ffdf(file="allPeptides.txt",
                                                  sep='\t', VERBOSE=FALSE,
                                                  colClasses=tempclasses,
                                                  next.rows=100000)

            allpeptides <- base::as.data.frame(allpeptides_save)
            allpeptides <- allpeptides[order(allpeptides$Mass), ]

            # Free some memory
            rm(allpeptides_save)
            gc()

            if(!is.null(sample_list)){
                # Check if all Raw.files are in sample_list and that all samples
                # in sample_list are in Raw.file
                all_samples_in_MaxQ <- sort(unique(allpeptides$Raw.file))

                # Any samples specified for IceR but not available in MaxQ data?
                if(any(sample_list %not in% all_samples_in_MaxQ)){
                    missing <- which(sample_list %not in% all_samples_in_MaxQ)
                    options(warn=1)
                    warning(paste0(paste(sample_list[missing], collapse = ","),
                            " is/are specified to be requantified by IceR but ",
                            "missing in MaxQuant results. This/these samples ",
                            "will be skipt by IceR."))
                    options(warn=-1)
                }

                # Any samples available in MaxQ but not specified for IceR?
                if(any(all_samples_in_MaxQ %not in% sample_list)){
                    missing <- which(all_samples_in_MaxQ %not in% sample_list)
                    options(warn=1)
                    warning(paste0(paste(all_samples_in_MaxQ[missing],
                                         collapse=",")," is/are available in ",
                                   "MaxQuant results but not specified to be ",
                                   "requantified by IceR. This/these samples ",
                                   "will be skipt by IceR."))
                    options(warn=-1)
                }
            }

            # Unify column names
            # RT
            # Possibly upper and lower case problem
            cnames <- colnames(allpeptides)
            if(length(which(cnames == "Retention.Length")) == 0){
                loc <- which(cnames == "Retention.length")

                if(length(loc) > 0){
                    colnames(allpeptides)[loc] <- "Retention.Length"
                }

                else{
                    allpeptides$Retention.Length <- 0
                }
            }

            if(MassSpec_mode == "Orbitrap"){
                # In Orbitrap allpeptides.txt the retention length is given in
                # sec instead of min like in e.g. evidence.txt
                aux <- allpeptides$Retention.Length / 60
                allpeptides$Retention.Length <- aux
            }

            # MSMS scan numbers
            cnames <- colnames(allpeptides)
            if(length(which(cnames == "MSMS.Scan.Numbers")) == 0){
                loc <- which(cnames == "MS.MS.scan.number")

                if(length(loc) > 0){
                    colnames(allpeptides)[loc] <- "MSMS.Scan.Numbers"
                }

                else{
                    allpeptides$MSMS.Scan.Numbers = 0
                }
            }

            # Some peptides are not available in allpeptides.txt thus we get
            # this additional information from the evidence.txt file
            evidence <- utils::read.csv(file="evidence.txt", sep='\t',
                                        header=TRUE)
            cnames <- colnames(evidence)

            if(any(cnames == "MS.MS.Scan.Number")){
                loc <- which(cnames == "MS.MS.Scan.Number")
                colnames(evidence)[loc] <- "MSMS.Scan.Numbers"
            }

            if(any(cnames == "MS.MS.scan.number")){
                loc <- which(cnames == "MS.MS.scan.number")
                colnames(evidence)[loc] <- "MSMS.Scan.Numbers"
            }

            if(any(cnames == "MSMS.Scan.Numbers")){
                loc <- which(cnames == "MSMS.Scan.Numbers")
                colnames(evidence)[loc] <- "MSMS.Scan.Numbers"
            }

            if(any(cnames == "Retention.length")){
                loc <- which(cnames == "Retention.length")
                colnames(evidence)[loc] <- "Retention.Length"
            }

            if(any(cnames == "Ion.mobility.length")){
                loc <- which(cnames == "Ion.mobility.length")
                colnames(evidence)[loc] <- "Ion.mobility.index.length"
            }

            if(any(cnames == "Number.of.scans")){
                loc <- which(cnames == "Number.of.scans")
                colnames(evidence)[loc] <- "Number.of.frames"
            }

            print(paste0(Sys.time(), " Read MaxQ results finished"))

            if(MassSpec_mode == "Orbitrap"){
                add_data <- base::as.data.frame(matrix(ncol=ncol(allpeptides),
                                                       nrow=nrow(evidence), NA))
                colnames(add_data) <- colnames(allpeptides)
                match_col_names <- match(colnames(allpeptides),
                                         colnames(evidence))
                jval <- as.integer(which(!is.na(match_col_names)))
                loc <- match_col_names[which(!is.na(match_col_names))]
                data.table::set(add_data, j=jval, value=evidence[, loc])

                # Now find rows (peptides + modification + sample) which are not
                # yet present in allpeptides.txt
                temp1 <- base::paste(add_data$Sequence, add_data$Modifications,
                                     add_data$Raw.file, sep="_")
                temp2 <- base::paste(allpeptides$Sequence,
                                     allpeptides$Modifications,
                                     allpeptides$Raw.file, sep="_")
                add_data <- add_data[which(temp1 %not in% temp2), ]
                allpeptides <- rbind(allpeptides, add_data)
                rm(temp1, temp2, add_data)
                gc()
            }

            if(MassSpec_mode == "TIMSToF"){
                add_data <- base::as.data.frame(matrix(ncol=ncol(allpeptides),
                                                       nrow=nrow(evidence), NA))
                colnames(add_data) <- colnames(allpeptides)
                match_col_names <- match(colnames(allpeptides),
                                         colnames(evidence))
                jval <- as.integer(which(!is.na(match_col_names)))
                loc <- match_col_names[which(!is.na(match_col_names))]
                data.table::set(add_data, j=jval, value=evidence[, loc])

                # In TIMS allpeptides.txt column sequence and modification is
                # missing thus take evidence table completeley
                add_data$Sequence <- evidence$Sequence
                add_data$Modifications <- evidence$Modifications
                add_data$Proteins <- evidence$Proteins
                add_data$Score <- evidence$Score

                allpeptides$Sequence <- ""
                allpeptides$Modifications <- ""
                allpeptides$Proteins <- ""
                allpeptides$Score <- NA
                allpeptides$Score <- as.numeric(allpeptides$Score)

                allpeptides <- rbind(allpeptides, add_data)
                rm(add_data)
                gc()

                # Convert Ion.mobility.index into 1/K0
                fit <- stats::lm(evidence$X1.K0~evidence$Ion.mobility.index)
                allpeptides$inverse_K0 <- ((allpeptides$Ion.mobility.index
                                            * fit$coefficients[2])
                                           + fit$coefficients[1])
                val <- allpeptides$Ion.mobility.index.length / 1000
                allpeptides$inverse_K0_length <- val
            }

            # If no specific sample list is defined which should be used, use
            # all raw files
            if(is.null(sample_list)){
                sample_list <- sort(unique(allpeptides$Raw.file))
            }

            # Clean allpeptides table to reduce required memory space
            if(MassSpec_mode == "Orbitrap"){
                if(any(colnames(allpeptides) == "Resolution")){
                    usecols <- c("Raw.file", "Charge", "m.z", "Mass",
                                 "Uncalibrated.m.z", "Resolution",
                                 "Max.intensity.m.z.0", "Retention.time",
                                 "Retention.Length", "MS.MS.IDs", "Sequence",
                                 "Modifications", "Proteins", "Score",
                                 "MSMS.Scan.Numbers", "Intensity")
                    allpeptides <- allpeptides[, usecols]
                }

                else{
                    usecols <- c("Raw.file", "Charge", "m.z", "Mass",
                                 "Uncalibrated.m.z", "Max.intensity.m.z.0",
                                 "Retention.time", "Retention.Length",
                                 "MS.MS.IDs", "Sequence", "Modifications",
                                 "Proteins", "Score", "MSMS.Scan.Numbers",
                                 "Intensity")
                    allpeptides <- allpeptides[, usecols]
                    allpeptides$Resolution <- 0
                    print("MS resolution seems to be missing in MaxQ outputs !")
                }
            }

            if(MassSpec_mode == "TIMSToF"){
                usecols <- c("Raw.file", "Charge", "m.z", "Mass",
                             "Retention.time", "Retention. Length", "Sequence",
                             "Modifications", "Proteins", "Score",
                             "MSMS.Scan.Numbers", "Intensity", "inverse_K0",
                             "inverse_K0_length")
                allpeptides <- allpeptides[, usecols]
                # No resolution information available
                allpeptides$Resolution <- 0
            }

            if(remove_contaminants == TRUE){
                exclude <- which(grepl("CON", allpeptides$Proteins)
                                 | grepl("REV", allpeptides$Proteins))
            }

            if(remove_contaminants == FALSE){
                exclude <- which(grepl("REV", allpeptides$Proteins))
            }

            if(length(exclude) > 0){
                # Remove potential reverse and contaminant peptides
                allpeptides <- allpeptides[-exclude, ]
            }

            # Add calibrated RT to all peptides
            usecols <- c("Sequence", "Raw.file", "Charge", "Modifications",
                         "Calibrated.retention.time")
            temp_evidence <- evidence[, usecols]
            rt <- temp_evidence$Calibrated.retention.time
            lst <- list(Sequence=evidence$Sequence, Raw.file=evidence$Raw.file,
                        Charge=evidence$Charge,
                        Modifications=evidence$Modifications)
            temp_evidence <- stats::aggregate(rt, lst, FUN=mean, na.rm=TRUE)
            colnames(temp_evidence)[5] <- "Calibrated.retention.time"
            usecols <- c("Sequence"="Sequence", "Raw.file"="Raw.file",
                         "Charge"="Charge", "Modifications"="Modifications")
            temp <- dplyr::left_join(allpeptides, temp_evidence, by=usecols)

            c1 <- is.na(temp$Calibrated.retention.time)
            c2 <- temp$Sequence != " "
            c3 <- temp$Sequence != ""
            missing_RT_calibration_indices <- which(c1 & c2 & c3)

            if(length(missing_RT_calibration_indices) > 0){
                for(ind in missing_RT_calibration_indices){
                    c1 <- temp$Sequence == temp$Sequence[ind]
                    c2 <- temp$Modifications == temp$Modifications[ind]
                    sub <- temp[which(c1 & c2), ]
                    aux <- stats::median(sub$Calibrated.retention.time,
                                         na.rm=TRUE)
                    temp$Calibrated.retention.time[ind] <- aux
                }
            }

            # If any peptide feature still doesn?t have a valid calibrated
            # RT --> just use observed RT
            c1 <- is.na(temp$Calibrated.retention.time)
            c2 <- temp$Sequence != " "
            c3 <- temp$Sequence != ""
            loc <- which(c1 & c2 & c3)
            temp$Calibrated.retention.time[loc] <- temp$Retention.time[loc]
            allpeptides <- temp
            print(paste0(Sys.time(), " Determine deviations of m/z from truth"))

            # Determine isotope corrected m/z
            if(MassSpec_mode == "Orbitrap"){
                # Calculate deviations between samples
                add <- base::as.data.frame(matrix(nrow=nrow(allpeptides),
                                                  ncol=3))
                colnames(add) <- c("isotope_at_max_int",
                                   "isotope_corrected_mz_at_max_int",
                                   "delta_mz_to_mz_at_max_int")
                add[, 1] <- as.numeric(add[, 1])
                add[, 2] <- as.numeric(add[, 2])
                add[, 3] <- as.numeric(add[, 3])

                isotopes <- matrix(nrow=nrow(allpeptides), ncol=4)
                isotopes[, 1] <- allpeptides$m.z
                isotopes[, 2] <- ((allpeptides$m.z * allpeptides$Charge)
                                  + (1 * 1.002054)) / allpeptides$Charge
                isotopes[, 3] <- ((allpeptides$m.z * allpeptides$Charge)
                                  + (2 * 1.002054)) / allpeptides$Charge
                isotopes[, 4] <- ((allpeptides$m.z * allpeptides$Charge)
                                  + (3 * 1.002054)) / allpeptides$Charge

                mx <- nrow(allpeptides)
                title <- "Determine deviations of true m/z from observed m/z"
                label <- base::paste(round(0 / max * 100, 0), "% done")
                pb <- tcltk::tkProgressBar(title=title, label=label, min=0,
                                           max=mx, width=300)
                start_time <- Sys.time()
                updatecounter <- 0
                time_require <- 0

                for(i in 1:nrow(allpeptides)){
                    if(!is.na(allpeptides$Max.intensity.m.z.0[i])){
                        mz_max_int <- allpeptides$Max.intensity.m.z.0[i]
                        deltas <- abs(isotopes[i, ] - mz_max_int)
                        mn <- min(deltas, na.rm=TRUE)
                        closest_iso <- which(deltas == mn) - 1
                        # +1 or +2 or +3 isotope shows highest intensity
                        # Correct back to isotope +0 m/z
                        if(closest_iso > 0){
                            iso_mult <- closest_iso * 1.002054
                            mz_max_int <- ((mz_max_int * allpeptides$Charge[i])
                                           - (iso_mult)) / allpeptides$Charge[i]
                        }

                        delta_mz_to_mz_at_max_int <- (mz_max_int
                                                      - allpeptides$m.z[i])
                        val <- c(closest_iso, mz_max_int,
                                 delta_mz_to_mz_at_max_int)
                        data.table::set(add, as.integer(i), as.integer(1:3),
                                        value=as.list(val))
                    }

                    updatecounter <- updatecounter + 1

                    if(updatecounter >= 100){
                        time_elapsed <- difftime(Sys.time(),
                                                 start_time,units="secs")
                        time_require <- ((time_elapsed / (i / mx))
                                         * (1 - (i / mx)))
                        td <- lubridate::seconds_to_period(time_require)
                        time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                                lubridate::minute(td),
                                                round(lubridate::second(td),
                                                      digits=0))
                        updatecounter <- 0
                        label <- base::paste(round(i / mx * 100, 0),
                                             " % done (", i, "/", max,
                                             ", Time require: ", time_require,
                                             ")", sep="")
                        tcltk::setTkProgressBar(pb, i, label=label)
                    }
                }

                close(pb)
                allpeptides <- cbind(allpeptides, add)

                # Evidence table
                add <- base::as.data.frame(matrix(nrow=nrow(evidence), ncol=3))
                colnames(add) <- c("isotope_at_max_int",
                                   "isotope_corrected_mz_at_max_int",
                                   "delta_mz_to_mz_at_max_int")
                add[, 1] <- as.numeric(add[, 1])
                add[, 2] <- as.numeric(add[, 2])
                add[, 3] <- as.numeric(add[, 3])

                isotopes <- matrix(nrow=nrow(evidence), ncol=4)
                isotopes[, 1] <- evidence$m.z
                isotopes[, 2] <- ((evidence$m.z * evidence$Charge)
                                  + (1 * 1.002054)) / evidence$Charge
                isotopes[, 3] <- ((evidence$m.z * evidence$Charge)
                                  + (2 * 1.002054)) / evidence$Charge
                isotopes[, 4] <- ((evidence$m.z * evidence$Charge)
                                  + (3 * 1.002054)) / evidence$Charge

                mx <- nrow(evidence)
                title <- "Determine deviations of true m/z from observed m/z"
                label <- base::paste(round(0 / mx * 100, 0), "% done")
                pb <- tcltk::tkProgressBar(title=title, label=label, min=0,
                                           max=mx, width=300)
                start_time <- Sys.time()
                updatecounter <- 0
                time_require <- 0

                for(i in 1:nrow(evidence)){
                    if(!is.na(evidence$Max.intensity.m.z.0[i])){
                        mz_max_int <- evidence$Max.intensity.m.z.0[i]
                        deltas <- abs(isotopes[i, ] - mz_max_int)
                        mn <- min(deltas, na.rm=TRUE)
                        closest_iso <- which(deltas == mn) - 1

                        # +1 or +2 or +3 isotope shows highest intensity
                        # Correct back to isotope +0 m/z
                        if(closest_iso > 0){
                            iso_mult <- closest_iso * 1.002054
                            mz_max_int <- ((mz_max_int * evidence$Charge[i])
                                           - (iso_mult)) / evidence$Charge[i]
                        }

                        delta_mz_to_mz_at_max_int <- (mz_max_int
                                                      - evidence$m.z[i])
                        val <- c(closest_iso, mz_max_int,
                                 delta_mz_to_mz_at_max_int)
                        data.table::set(add, as.integer(i), as.integer(1:3),
                                        value=as.list(val))
                    }

                    updatecounter <- updatecounter + 1

                    if(updatecounter >= 100){
                        time_elapsed <- difftime(Sys.time(), start_time,
                                                 units="secs")
                        time_require <- ((time_elapsed / (i / mx))
                                         * (1 - (i / mx)))
                        td <- lubridate::seconds_to_period(time_require)
                        time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                                lubridate::minute(td),
                                                round(lubridate::second(td),
                                                      digits=0))
                        updatecounter <- 0
                        label <- base::paste(round(i / mx * 100, 0),
                                             " % done (", i, "/", mx,
                                             ", Time require: ", time_require,
                                             ")", sep="")
                        tcltk::setTkProgressBar(pb, i, label=label)
                    }
                }

                close(pb)
                evidence <- cbind(evidence, add)
            }

            # Finally save prepared data
            setwd(path_to_output)
            save(allpeptides, file="Temporary_files/allPeptides.RData")
            save(evidence, file="Temporary_files/evidence.RData")
        }

        # Check up-front that all samples were measured with comparable gradient
        # lengths
        max_RTs <- stats::aggregate(allpeptides$Retention.time,
                                    by=list(Raw.file=allpeptides$Raw.file),
                                    mx, na.rm=TRUE)
        max_RTs_range <- c(round(min(max_RTs$x, na.rm=TRUE), digits=0),
                           round(max(max_RTs$x, na.rm=TRUE), digits=0))

        if(max_RTs_range[2] - max_RTs_range[1] > 1){
            stop(paste0("Gradient lengths between samples differ too much!!! ",
                        "Minimal observed gradient length: ", max_RTs_range[1],
                        ", Maximal observed gradient length: ",
                        max_RTs_range[2], ". IceR can only reliably process ",
                        "data acquired with same gradient lengths!!!"))
        }

        # Here relevant qc data is stored and finally saved as RData which can
        # be used for re-generating plots
        QC_data <- list()

        grDevices::pdf(base::paste("Temporary_files/QC_plots",
                                   output_file_names_add, ".pdf", sep=""))

        RT_calibration <- TRUE
        mz_calibration <- TRUE

        if(MassSpec_mode == "TIMSToF"){
            IM_calibration <- TRUE
        }

        else{
            IM_calibration <- FALSE
        }

        if(MassSpec_mode == "Orbitrap"){
            usecols <- c("Retention.time.calibration",
                         "Uncalibrated...Calibrated.m.z..Da.")
            QC_data[["MaxQ_calibrations"]] <- evidence[, usecols]
        }

        allpeptides <- allpeptides[order(allpeptides$m.z), ]
        rownames(allpeptides) <- c(1:nrow(allpeptides))

        c1 <- Sequence != " "
        c2 <- Sequence != ""
        c3 <- Modifications == "Unmodified"

        if(only_unmodified_peptides == TRUE){
            identified_ions <- subset(allpeptides, c1 & c2 & c3)
        }

        else{
            identified_ions <- subset(allpeptides, c1 & c2)
        }

        identified_ions <- identified_ions[order(identified_ions$Sequence,
                                                 identified_ions$Intensity,
                                                 decreasing=TRUE), ]

        identified_ions$Sequence <- as.character(identified_ions$Sequence)
        unique_peptides <- unique(identified_ions$Sequence)

        if(MassSpec_mode == "Orbitrap"){
            windows <- base::as.data.frame(matrix(ncol=11,
                                                  nrow=nrow(identified_ions)))
            colnames(windows) <- c("Peptide", "Mod", "Charge", "mean_m.z",
                                   "sd_m.z", "mean_RT", "sd_RT", "num_ions",
                                   "mean_RT_calibrated", "sd_RT_calibrated",
                                   "sd_m.z_uncalibrated")

            windows$Peptide <- as.character(windows$Peptide)
            windows$Mod <- as.character(windows$Mod)
            windows$Charge <- as.numeric(windows$Charge)
            windows$mean_m.z <- as.numeric(windows$mean_m.z)
            windows$sd_m.z <- as.numeric(windows$sd_m.z)
            windows$mean_RT <- as.numeric(windows$mean_RT)
            windows$sd_RT <- as.numeric(windows$sd_RT)
            windows$num_ions <- as.numeric(windows$num_ions)
            windows$mean_RT_calibrated <- as.numeric(windows$mean_RT_calibrated)
            windows$sd_RT_calibrated <- as.numeric(windows$sd_RT_calibrated)
            aux <- as.numeric(windows$sd_m.z_uncalibrated)
            windows$sd_m.z_uncalibrated <- aux
        }

        if(MassSpec_mode == "TIMSToF"){
            windows <- base::as.data.frame(matrix(ncol=13,
                                                  nrow=nrow(identified_ions)))
            colnames(windows) <- c("Peptide", "Mod", "Charge", "mean_m.z",
                                   "sd_m.z", "mean_RT", "sd_RT", "num_ions",
                                   "mean_RT_calibrated", "sd_RT_calibrated",
                                   "sd_m.z_uncalibrated", "mean_IM", "sd_IM")

            windows$Peptide <- as.character(windows$Peptide)
            windows$Mod <- as.character(windows$Mod)
            windows$Charge <- as.numeric(windows$Charge)
            windows$mean_m.z <- as.numeric(windows$mean_m.z)
            windows$sd_m.z <- as.numeric(windows$sd_m.z)
            windows$mean_RT <- as.numeric(windows$mean_RT)
            windows$sd_RT <- as.numeric(windows$sd_RT)
            windows$num_ions <- as.numeric(windows$num_ions)
            windows$mean_RT_calibrated <- as.numeric(windows$mean_RT_calibrated)
            windows$sd_RT_calibrated <- as.numeric(windows$sd_RT_calibrated)
            aux <- as.numeric(windows$sd_m.z_uncalibrated)
            windows$sd_m.z_uncalibrated <- aux
            windows$mean_IM <- as.numeric(windows$mean_IM)
            windows$sd_IM <- as.numeric(windows$sd_IM)
        }

        # Remove outlier samples where RT of a peptide is very different from
        # all other samples
        # Deviation should not be larger than 5 % of the total chromatographic
        # retention length
        outlier_RT_deviation <- (max(allpeptides$Retention.time, na.rm=TRUE)
                                 / 100) * 5
        print(paste0(Sys.time(), " Determine matching parameter windows"))

        mx <- length(unique_peptides)
        label <- base::paste(round(0 / mx * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title="Determine matching parameter windows",
                                   label=label, min=0, max=mx, width=300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0
        start <- 1
        max_ind <- nrow(identified_ions)
        count_features <- 0

        for(i in 1:length(unique_peptides)){
            ind <- start

            while(TRUE){
                ind <- ind + 100

                if(ind < max_ind){
                    if(identified_ions[ind, "Sequence"] != unique_peptides[i]){
                        break
                    }
                }

                else{
                    ind <- max_ind
                    break
                }
            }

            sub <- identified_ions[start:ind, ]

            # Reset for next start index
            bins <- (ind - start) / 100
            start <- ifelse(bins == 1, start, start + ((bins - 1) * 100))

            sub <- sub[which(sub$Sequence == unique_peptides[i]), ]

            if(length(which(!is.na(sub$Intensity))) > 0){
                sub<-sub[which(!is.na(sub$Intensity)), ]
            }

            sub <- sub[!duplicated(sub$Raw.file), ]

            for(c in unique(sub$Charge)){
                sub1 <- sub[which(sub$Charge == c), ]

                for(m in unique(sub1$Modifications)){
                    sub1 <- sub[which(sub$Charge == c
                                      & sub$Modifications == m), ]
                    count_features <- count_features + 1

                    median_RT <- stats::median(sub1$Retention.time, na.rm=TRUE)
                    aux <- median_RT + outlier_RT_deviation
                    remove <- which(sub1$Retention.time > aux)

                    if(length(remove) > 0){
                        sub1 <- sub1[-remove, ]
                    }

                    # Determine mean m/z at peak maximum and mean RT
                    if(MassSpec_mode == "Orbitrap"){
                        lst <- c(c, mean(sub1$m.z, na.rm=TRUE),
                                 stats::sd(sub1$m.z, na.rm=TRUE),
                                 mean(sub1$Retention.time, na.rm=TRUE),
                                 stats::sd(sub1$Retention.time, na.rm=TRUE),
                                 nrow(sub1),
                                 mean(sub1$Calibrated.retention.time,
                                      na.rm=TRUE),
                                 stats::sd(sub1$Calibrated.retention.time,
                                           na.rm=TRUE),
                                 stats::sd(sub1$isotope_corrected_mz_at_max_int,
                                           na.rm=TRUE))
                        data.table::set(windows, as.integer(count_features),
                                        as.integer(c(3:11)), value=as.list(lst))
                    }

                    if(MassSpec_mode == "TIMSToF"){
                        lst <- c(c, mean(sub1$m.z, na.rm=TRUE),
                                 stats::sd(sub1$m.z, na.rm=TRUE),
                                 mean(sub1$Retention.time, na.rm=TRUE),
                                 stats::sd(sub1$Retention.time, na.rm=TRUE),
                                 nrow(sub1),
                                 mean(sub1$Calibrated.retention.time,
                                      na.rm=TRUE),
                                 stats::sd(sub1$Calibrated.retention.time,
                                           na.rm=TRUE),
                                 stats::sd(sub1$isotope_corrected_mz_at_max_int,
                                           na.rm=TRUE),
                                 mean(sub1$inverse_K0, na.rm=TRUE),
                                 stats::sd(sub1$inverse_K0, na.rm=TRUE))
                        data.table::set(windows, as.integer(count_features),
                                        as.integer(c(3:13)), value=as.list(lst))
                    }

                    lst <- as.list(c(unique_peptides[i], m))
                    data.table::set(windows, as.integer(count_features),
                                    as.integer(1:2), value=lst)
                }
            }

            updatecounter <- updatecounter + 1

            if(updatecounter >= 10){
                time_elapsed <- difftime(Sys.time(), start_time, units="secs")
                time_require <- (time_elapsed / (i / mx)) * (1 - (i / mx))
                td <- lubridate::seconds_to_period(time_require)
                time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                        lubridate::minute(td),
                                        round(lubridate::second(td), digits=0))

                updatecounter <- 0
                label <- base::paste(round(i / mx * 100, 0), " % done (", i,
                                     "/", mx, ", Time require: ", time_require,
                                     ")", sep="")
                tcltk::setTkProgressBar(pb, i, label=label)
            }
        }

        close(pb)

        windows <- windows[1:count_features, ]

        # Generate some QC plots indicating variation of RT and m/z (and IM) for
        # identified features across samples
        graphics::boxplot(windows$sd_m.z_uncalibrated,
                          main="Standard deviation of peptide features m/z",
                          outline=FALSE, ylab="StDev of m/z [Da]")
        graphics::boxplot(windows$sd_RT,
                          main="Standard deviation of peptide features RT",
                          outline=FALSE, ylab="StDev of RT [min]")

        # Return some values to the console
        rt_dev <- median(windows$sd_RT, na.rm=TRUE)
        mz_dev <- median(windows$sd_m.z_uncalibrated, na.rm=TRUE)
        print(paste0("Median mz-deviation for peptides between samples: ",
                     round(mz_dev, digits=6), " Da"))
        print(paste0("Median RT-deviation for peptides between samples: ",
                     round(rt_dev, digits=2), " min"))

        options(warn=1)

        if(rt_dev > 1){
            warning(pate0("Median RT-deviation between samples is larger than ",
                          "1 min which could indicate that chromatography is ",
                          "not stable. IceR might not be able to process the ",
                          "data correctly!!!"))
        }

        if(mz_dev > 0.01){
            warning(paste0("Median m/z-deviation between samples is larger ",
                           "than 0.01 Da which could indicate that Mass ",
                           "analyzer is not stable. IceR might not be able to ",
                           "process the data correctly!!!"))
        }

        options(warn=-1)

        if(MassSpec_mode == "TIMSToF"){
            title <- "Standard deviation of peptide features inverse K0"
            graphics::boxplot(windows$sd_IM, main=title, outline=FALSE)
        }

        if(is.na(mz_window)){
            # Based on boxplots
            aux <- grDevices::boxplot.stats(windows$sd_m.z)$`stats`[5]

            if(!is.na(min_mz_window)){
                if(aux > min_mz_window){
                    borders_m.z <- c(-aux, aux) # Upper whisker
                }

                else{
                    borders_m.z <- c(-min_mz_window, min_mz_window)
                }
            }

            else{
                borders_m.z <- c(-aux, aux) # Upper whisker
            }
        }

        else{
            # Used defined parameter
            borders_m.z <- c(-mz_window, mz_window)
        }

        if(is.na(RT_window)){
            # Based on boxplots
            aux <- grDevices::boxplot.stats(windows$sd_RT)$`stats`[5]

            if(!is.na(min_RT_window)){
                if(aux > min_RT_window){
                    borders_RT <- c(-aux, aux) # 75% quantil based on boxplots
                }

                else{
                    borders_RT <- c(-min_RT_window, min_RT_window)
                }
            }

            else{
                borders_RT <- c(-aux, aux) # 75% quantil based on boxplots
            }
        }

        else{
            # Use defined parameter
            borders_RT <- c(-RT_window, RT_window)
        }

        if(MassSpec_mode == "TIMSToF"){
            if(is.na(IM_window)){
                # Based on boxplots
                aux <- grDevices::boxplot.stats(windows$sd_IM)$`stats`[5]

                if(!is.na(min_IM_window)){
                    if(aux > min_IM_window){
                        # 75% quantil based on boxplots
                        borders_IM <- c(-aux, aux)
                    }

                    else{
                        borders_IM <- c(-min_IM_window, min_IM_window)
                    }
                }

                else{
                    borders_IM <- c(-aux, aux) # 75% quantil based on boxplots
                }
            }

            else{
                # Use defined parameter
                borders_IM <- c(-IM_window, IM_window)
            }
        }

        borders_RT_use_save =<-borders_RT

        if(MassSpec_mode == "Orbitrap"){
            lst <- list(RT_window=borders_RT, mz_window=borders_m.z)
            QC_data[["Feature_alignment_windows"]] <- lst
        }

        if(MassSpec_mode == "TIMSToF"){
            lst <- list(RT_window=borders_RT, mz_window=borders_m.z,
                        IM_window=borders_IM)
            QC_data[["Feature_alignment_windows"]] <- lst
        }


        # Determine for how many windows we expect an overlap of ions for
        # different peptide sequences based on the chosen parameters
        windows <- windows[order(windows$mean_m.z), ]
        max_i <- nrow(windows)
        windows$overlap <- 0

        mx <- nrow(windows)
        title <-"Determine overlapping peptide features"
        label <- base::paste(round(0 / mx * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title=title, label=label, min=0, max=mx,
                                   width=300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0

        if(MassSpec_mode == "Orbitrap"){
            for(i in 1:nrow(windows)){
                start <- ifelse((i - 5) < 1, 1, (i - 5))
                end <- ifelse((i + 5) > max_i, max_i, (i + 5))
                sub <- windows[start:end, ]
                c1 <- sub$mean_m.z >= windows$mean_m.z[i] + borders_m.z[1]
                c2 <- sub$mean_m.z <= windows$mean_m.z[i] + borders_m.z[2]
                c3 <- sub$mean_RT >= windows$mean_RT[i] + borders_RT[1]
                c4 <- sub$mean_RT <= windows$mean_RT[i] + borders_RT[2]
                c5 <- sub$Charge == windows$Charge[i]
                sub <- sub[which(c1 & c2 & c3 & c4 & c5), ]

                data.table::set(windows, as.integer(i), 12L,
                                value=(nrow(sub) - 1))

                updatecounter <- updatecounter + 1

                if(updatecounter >= 100){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / mx)) * (1 - (i / mx))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))
                    updatecounter <- 0
                    label <- base::paste(round(i / mx * 100, 0), " % done (", i,
                                         "/", mx, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }
        }

        if(MassSpec_mode == "TIMSToF"){
            for(i in 1:nrow(windows)){
                start <- ifelse((i - 5) < 1, 1, (i - 5))
                end <- ifelse((i + 5) > max_i, max_i, (i + 5))
                sub <- windows[start:end, ]
                c1 <- sub$mean_m.z >= windows$mean_m.z[i] + borders_m.z[1]
                c2 <- sub$mean_m.z <= windows$mean_m.z[i] + borders_m.z[2]
                c3 <- sub$mean_RT >= windows$mean_RT[i] + borders_RT[1]
                c4 <- sub$mean_RT <= windows$mean_RT[i] + borders_RT[2]
                c5 <- sub$mean_IM >= windows$mean_IM[i] + borders_IM[1]
                c6 <- sub$mean_IM <= windows$mean_IM[i] + borders_IM[2]
                c7 <- sub$Charge == windows$Charge[i]
                sub <- sub[which(c1 & c2 & c3 & c4 & c5 & c6 & c7), ]
                data.table::set(windows, as.integer(i), 12L,
                                value=(nrow(sub) - 1))
                updatecounter <- updatecounter + 1

                if(updatecounter >= 100){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / mx)) * (1 - (i / mx))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))
                    updatecounter <- 0
                    label <- base::paste(round(i / mx * 100, 0), " % done (", i,
                                         "/", mx, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }
        }

        close(pb)

        p <- Barplots(plyr::count(windows$overlap)[, 2],
                      main="Number of expected feature overlaps", ylab="Count",
                      Name=plyr::count(windows$overlap)[, 1], xlab="Overlaps",
                      AvgLine=FALSE)
        QC_data[["Alignment_deviations_overlap"]] <- windows

        # Perform alignment of peptide features
        allpeptides <- allpeptides[order(allpeptides$m.z), ]
        rownames(allpeptides) <- c(1:nrow(allpeptides))
        allpeptides$Sequence <- as.character(allpeptides$Sequence)
        allpeptides$Modifications <- as.character(allpeptides$Modifications)
        allpeptides$Proteins <- base::gsub(";", "|", allpeptides$Proteins)
        allpeptides$Raw.file <- as.character(allpeptides$Raw.file)

        # Presubset all ions based on charge, here use only rows of all ions
        # which are coming from samples which should be also used
        # Also create a library of ions per charge state per mz_window of 0.5 Da
        print(paste0(Sys.time(), " Prepare for peptide feature matching - ",
                     "Indexing all ions"))
        rownames(allpeptides) <- c(1:nrow(allpeptides))
        allpeptides_frag <- list()
        allpeptides_frag_indices_per_mz_window <- list()

        for(i in 1:max(allpeptides$Charge)){
            loc <- which(allpeptides$Charge == i)

            if(length(loc) > 0){
                allpeptides_frag[[i]] <- allpeptides[loc, ]
                min_max_mz <- c(floor(min(allpeptides_frag[[i]]$m.z)),
                                ceiling(max(allpeptides_frag[[i]]$m.z)))
                indices <- matrix(ncol=4, nrow=min_max_mz[2] - min_max_mz[1])
                indices[,1] <- min_max_mz[1]:(min_max_mz[2] - 1)
                indices[,2] <- (min_max_mz[1] + 1):(min_max_mz[2])

                for(j in 1:nrow(indices)){
                    indx <- which(allpeptides_frag[[i]]$m.z >= indices[j, 1]
                                  & allpeptides_frag[[i]]$m.z <= indices[j, 2])
                    if(length(indx) > 0){
                        indices[j, 3:4] <- c(min(indx), max(indx))
                    }

                    else{
                        indices[j, 3:4] <- c(NA, NA)
                    }
                }

                allpeptides_frag_indices_per_mz_window[[i]] <- indices
            }

            else{
                allpeptides_frag[[i]] <- NA
                allpeptides_frag_indices_per_mz_window[[i]] <- NA
            }
        }

        rm(indices)
        gc()

        # Here already used ions will be marked
        temp_data <- base::as.data.frame(matrix(ncol=1, nrow=nrow(allpeptides)))
        temp_data[, 1] <- as.numeric(temp_data[, 1])
        temp <- subset(allpeptides, Sequence != " " & Sequence != "")

        if(MassSpec_mode == "Orbitrap"){
            nrow <- length(unique(base::paste(temp$Sequence, temp$Charge,
                                              temp$Modifications)))
            features <- base::as.data.frame(matrix(ncol=21, nrow=nrow))
            colnames(features) <- c("Feature_name", "m.z", "RT", "Sequence",
                                    "Protein", "MSMS.Scan.Numbers",
                                    "m.z_range_min", "m.z_range_max",
                                    "RT_range_min", "RT_range_max",
                                    "ion_indices", "count_ion_indices", "Mass",
                                    "Charge", "mean_Scores", "num_matches",
                                    "Modifications", "RT_length", "Observed_RT",
                                    "Observed_mz", "Observed_score")

            features <- data.table::as.data.table(features)
            features$Feature_name <- as.character(features$Feature_name)
            features$Mass <- as.numeric(features$Mass)
            features$RT <- as.numeric(features$RT)
            features$Sequence <- as.character(features$Sequence)
            features$Protein <- as.character(features$Protein)
            scnum <- features$MSMS.Scan.Numbers
            features$MSMS.Scan.Numbers <- as.character(scnum)
            features$m.z_range_min <- as.numeric(features$m.z_range_min)
            features$m.z_range_max <- as.numeric(features$m.z_range_max)
            features$RT_range_min <- as.numeric(features$RT_range_min)
            features$RT_range_max <- as.numeric(features$RT_range_max)
            features$ion_indices <- as.character(features$ion_indices)
            features$count_ion_indices <- as.numeric(features$count_ion_indices)
            features$m.z <- as.numeric(features$m.z)
            features$Charge <- as.numeric(features$Charge)
            features$mean_Scores <- as.character(features$mean_Scores)
            features$num_matches <- as.character(features$num_matches)
            features$Modifications <- as.character(features$Modifications)
            features$RT_length <- as.numeric(features$RT_length)
            features$Observed_RT <- as.character(features$Observed_RT)
            features$Observed_mz <- as.character(features$Observed_mz)
            features$Observed_score <- as.character(features$Observed_score)
        }

        if(MassSpec_mode == "TIMSToF"){
            nrow <- length(unique(base::paste(temp$Sequence, temp$Charge,
                                              temp$Modifications)))
            features <- base::as.data.frame(matrix(ncol=26, nrow=nrow))
            colnames(features) <- c("Feature_name", "m.z", "RT", "Sequence",
                                    "Protein", "MSMS.Scan.Numbers",
                                    "m.z_range_min", "m.z_range_max",
                                    "RT_range_min", "RT_range_max",
                                    "ion_indices", "count_ion_indices", "Mass",
                                    "Charge", "mean_Scores", "num_matches",
                                    "Modifications", "RT_length", "Observed_RT",
                                    "Observed_mz", "Observed_score", "Inv_K0",
                                    "Inv_K0_length", "Inv_K0_range_min",
                                    "Inv_K0_range_max", "Observed_IM")

            features <- data.table::as.data.table(features)
            features$Feature_name <- as.character(features$Feature_name)
            features$Mass <- as.numeric(features$Mass)
            features$RT <- as.numeric(features$RT)
            features$Sequence <- as.character(features$Sequence)
            features$Protein <- as.character(features$Protein)
            scnum <- features$MSMS.Scan.Numbers
            features$MSMS.Scan.Numbers <- as.character(scnum)
            features$m.z_range_min <- as.numeric(features$m.z_range_min)
            features$m.z_range_max <- as.numeric(features$m.z_range_max)
            features$RT_range_min <- as.numeric(features$RT_range_min)
            features$RT_range_max <- as.numeric(features$RT_range_max)
            features$ion_indices <- as.character(features$ion_indices)
            features$count_ion_indices <- as.numeric(features$count_ion_indices)
            features$m.z <- as.numeric(features$m.z)
            features$Charge <- as.numeric(features$Charge)
            features$mean_Scores <- as.character(features$mean_Scores)
            features$num_matches <- as.character(features$num_matches)
            features$Modifications <- as.character(features$Modifications)
            features$RT_length <- as.numeric(features$RT_length)
            features$Observed_RT <- as.character(features$Observed_RT)
            features$Observed_mz <- as.character(features$Observed_mz)
            features$Observed_score <- as.character(features$Observed_score)
            features$Inv_K0 <- as.numeric(features$Inv_K0)
            features$Inv_K0_length <- as.numeric(features$Inv_K0_length)
            features$Inv_K0_range_min <- as.numeric(features$Inv_K0_range_min)
            features$Inv_K0_range_max <- as.numeric(features$Inv_K0_range_max)
            features$Observed_IM <- as.character(features$Observed_IM)
        }

        # Get ordering of allpeptides based on sequence
        order_by_Sequence <- order(allpeptides$Sequence, decreasing=TRUE)
        unique_peptides <- sort(unique_peptides, decreasing=TRUE)
        aux <- allpeptides$Sequence[which(allpeptides$Sequence != " "
                                          & allpeptides$Sequence != "")]
        num_rows_x_pept_seq <- plyr::count(aux)
        num_rows_x_pept_seq <- num_rows_x_pept_seq[order(num_rows_x_pept_seq$x,
                                                         decreasing=TRUE), ]

        max <- length(unique_peptides)
        label <- base::paste(round(0 / max * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title="Matching known features",
                                   label=label, min=0, max=max, width=300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0
        # Required for subsetting allpeptides per peptide sequence
        last_index <- 1
        count_features <- 0

        if(MassSpec_mode == "Orbitrap"){
            for(i in 1:length(unique_peptides)){
                loc <- last_index:(last_index - 1 + num_rows_x_pept_seq[i, 2])
                sub_peptide <- allpeptides[order_by_Sequence[loc], ]
                last_index <- last_index + num_rows_x_pept_seq[i, 2]
                loc <- which(sub_peptide$Sequence == unique_peptides[i])
                sub_peptide <- sub_peptide[loc, ]

                for(c in unique(sub_peptide$Charge)){
                    for(m in unique(sub_peptide$Modifications)){
                        c1 <- sub_peptide$Charge == c
                        c2 <- sub_peptide$Modifications == m
                        sub <- sub_peptide[which(c1 & c2), ]

                        if(nrow(sub) > 0){
                            sub <- sub[order(sub$Intensity,decreasing=TRUE), ]

                            # Only use features close to median feature over all
                            # samples for this peptide
                            median_m.z <- stats::median(sub$m.z, na.rm=TRUE)
                            median_RT <- stats::median(sub$Retention.time,
                                                       na.rm=TRUE)
                            d1 <- abs(sub$Retention.time - median_RT)
                            d2 <- abs(sub$m.z - median_m.z)
                            selection <- which(d1 < borders_RT[2]
                                               & d2 < borders_m.z[2])
                            sub <- sub[selection, ]

                            if(nrow(sub) > 0){
                                # Check if peptide was sequenced several times.
                                # If this is the case, use feature closest to
                                # median
                                if(any(duplicated(sub$Raw.file))){
                                    sub <- sub[!duplicated(sub$Raw.file), ]
                                }

                                median_m.z <- stats::median(sub$m.z, na.rm=TRUE)
                                median_RT <- stats::median(sub$Retention.time,
                                                           na.rm=TRUE)
                                mx <- max(sub$Retention.time, na.rm=TRUE)
                                mn <- min(sub$Retention.time, na.rm=TRUE)
                                RT_max_min <- mx - mn
                                too_large_RT_variation <- FALSE

                                # Get subset of allpeptides based on charge and
                                # mz window indexing
                                p <- allpeptides_frag_indices_per_mz_window[[c]]
                                mz1 <- median_m.z + borders_m.z[1]
                                mz2 <- median_m.z + borders_m.z[2]
                                start_ind <- which(p[, 1] > mz1)[1] - 1
                                end_ind <- which(p[, 2] > mz2)[1] + 1

                                if(is.na(end_ind) | end_ind > nrow(p)){
                                    end_ind <- nrow(p)
                                }

                                c1 <- is.na(p[start_ind, 3])
                                c2 <- !is.na(p[start_ind + 1, 3])

                                if(c1 & c2){
                                    start_ind <- p[start_ind + 1, 3]
                                    end_ind_temp <- p[end_ind, 4]

                                    if(is.na(end_ind_temp)){
                                        aux <- p[(end_ind + 1):nrow(p), 3]
                                        end_ind <- na.omit(aux)[1]
                                    }

                                    else{
                                        end_ind <- end_ind_temp
                                    }
                                }

                                c1 <- !is.na(p[start_ind, 3])
                                c2 <- is.na(p[end_ind, 4])

                                else if(c1 & c2){
                                    end_ind <- p[start_ind, 4]
                                    start_ind <- p[start_ind, 3]
                                }

                                c1 <- is.na(p[start_ind, 3])
                                c2 <- !is.na(p[end_ind, 3])

                                else if(c1 & c2){
                                    start_ind <- p[end_ind, 3]
                                    end_ind <- p[end_ind, 4]
                                }

                                else{
                                    start_ind <- p[start_ind, 3]
                                    end_ind <- p[end_ind, 4]
                                }

                                loc <- start_ind:end_ind
                                temp_sub <- allpeptides_frag[[c]][loc, ]
                                tmz <- temp_sub$m.z
                                c1 <- tmz >= median_m.z + borders_m.z[1]
                                c2 <- tmz <= median_m.z + borders_m.z[2]
                                trt <- temp_sub$Retention.time
                                c3 <- trt >= median_RT + borders_RT[1]
                                c4 <- trt <= median_RT + borders_RT[2]
                                c5 <- temp_sub$Sequence == unique_peptides[i]
                                c6 <- temp_sub$Charge == c
                                c7 <- temp_sub$Modifications == m
                                loc <- which(c1 & c2 & c3 & c4 | c5 & c6 & c7)
                                sub <- temp_sub[loc, ]
                                sub <- sub[order(sub$Intensity, sub$Sequence,
                                                 decreasing=TRUE), ]

                                # Check if selected MaxQ features in expected RT
                                # and mz window contains outliers (too high
                                # deviation in RT or mz compared to all other
                                # features)
                                box_stats_rt <- grDevices::boxplot.stats(
                                    sub$Retention.time
                                )$stats

                                box_stats_mz <- grDevices::boxplot.stats(
                                    sub$m.z
                                )$stats

                                c1 <- sub$Retention.time < box_stats_rt[2]
                                c2 <- sub$Retention.time > box_stats_rt[4]
                                c3 <- sub$m.z < box_stats_mz[2]
                                c4 <- sub$m.z > box_stats_mz[4]
                                outlier <- which(c1 | c2 | c3 | c4)

                                # Features with sequence identification will be
                                # never regarded as outlier
                                c <- sub$Sequence[outlier] == unique_peptides[i]

                                if(any(c)){
                                    outlier <- outlier[-which(c)]
                                }

                                # Any outlier detected? exclude them
                                if(length(outlier) > 0){
                                    sub <- sub[-outlier, ]
                                }

                                # Check if the fraction of features with not
                                # currently searched sequence is > 25 %
                                # If yes, keep all sequences and indicate this
                                # with multiple sequences for this feature
                                # otherwise remove other sequences
                                loc <- which(sub$Sequence != " "
                                             & sub$Sequence != "")
                                count_seqs <- plyr::count(sub$Sequence[loc])

                                loc <- which(count_seqs$x == unique_peptides[i])
                                aux <- 0.8 * sum(count_seqs$freq)

                                if(count_seqs$freq[loc] >= aux){
                                    c1 <- sub$Sequence == unique_peptides[i]
                                    c2 <- sub$Sequence == " "
                                    c3 <- sub$Sequence == ""
                                    sub <- sub[which(c1 | c2 | c3), ]
                                }

                                if(any(duplicated(sub$Raw.file))){
                                    sub <- sub[!duplicated(sub$Raw.file), ]
                                }

                                # Variation in RT for same peptide sequence
                                # between samples is too high
                                if(nrow(sub) == 0){
                                    borders_RT_save <- borders_RT
                                    too_large_RT_variation <- TRUE
                                    borders_RT <- c(-RT_max_min / 2,
                                                    RT_max_min / 2)
                                    tmz <- temp_sub$m.z
                                    c1 <- tmz >= median_m.z + borders_m.z[1]
                                    c2 <- tmz <= median_m.z + borders_m.z[2]
                                    trt <- temp_sub$Retention.time
                                    c3 <- trt >= median_RT + borders_RT[1]
                                    c4 <- trt <= median_RT + borders_RT[2]
                                    sub <- temp_sub[which(c1 & c2 & c3 & c4), ]

                                    # Still no match
                                    if(nrow(sub) == 0){
                                        too_large_RT_variation <- FALSE
                                        borders_RT <- borders_RT_save
                                    }
                                }

                                if(nrow(sub) > 0){
                                    count_features <- count_features + 1

                                    # New feature - give a new name
                                    name <- base::paste("Feature",
                                                        count_features, sep="_")

                                    nchr <- nchar(as.character(sub$Sequence))
                                    sequences_relevant <- which(nchr > 1)

                                    if(length(sequences_relevant) > 0){
                                        aux <- c(sub_peptide$Sequence[1],
                                                 sub[sequences_relevant,
                                                     "Sequence"])
                                        sequence <- base::paste(unique(aux),
                                                                collapse=";")
                                    }

                                    else{
                                        sequence <- sub_peptide$Sequence[1]
                                    }

                                    nchr <- nchar(as.character(sub$Proteins))
                                    proteins_relevant <- which(nchr > 1)

                                    if(length(proteins_relevant) > 0){
                                        aux <- c(sub_peptide$Proteins[1],
                                                 sub[proteins_relevant,
                                                     "Proteins"])
                                        protein <- base::paste(unique(aux),
                                                               collapse=";")
                                    }

                                    else{
                                        protein <- sub_peptide$Proteins[1]
                                    }

                                    chr <- as.character(sub$MSMS.Scan.Numbers)
                                    nchr <- nchar(chr) > 1
                                    msmsscansrelevant_relevant <- which(nchr)

                                    if(length(msmsscansrelevant_relevant) > 0){
                                        a <- sub[msmsscansrelevant_relevant,
                                                 "Raw.file"]
                                        b <- sub[msmsscansrelevant_relevant,
                                                 "MSMS.Scan.Numbers"]
                                        aux <- base::paste(a, "_msms_", b,
                                                           sep="")
                                        msmsscan <- base::paste(unique(aux),
                                                                collapse=";")
                                    }

                                    else{
                                        msmsscan <- ""
                                    }

                                    median_mass <- stats::median(sub$Mass,
                                                                 na.rm=TRUE)
                                    Charge <- c

                                    scores <- NULL
                                    num_matches <- NULL
                                    spsq <- stringr::str_split(sequence, ";")

                                    for(s in unique(unlist(spsq))){
                                        loc <- which(sub$Sequence == s)
                                        scores <- append(scores,
                                                         mean(sub$Score[loc],
                                                              na.rm=TRUE))
                                        num_matches <- append(num_matches,
                                                              nrow(sub[loc, ]))
                                    }

                                    Modifications <- ifelse(m != "Unmodified",
                                                            m, "")
                                    RT_length <- max(sub$Retention.Length,
                                                     na.rm=TRUE)
                                    lst <- list(Raw.file=sub$Raw.file)
                                    temp_RT <- stats::aggregate(
                                        sub$Retention.time,
                                        by=lst,
                                        FUN=mean,
                                        na.rm=TRUE
                                    )
                                    # Order all observed RT according to
                                    # sample_list
                                    Observed_RT <- base::paste(
                                        temp_RT[match(sample_list,
                                                      temp_RT$Raw.file), 2],
                                        collapse=";"
                                    )
                                    temp_mz <- stats::aggregate(
                                        sub$isotope_corrected_mz_at_max_int,
                                        by=lst,
                                        FUN=mean,
                                        na.rm=TRUE
                                    )
                                    # Order all observed mz according to
                                    # sample_list
                                    Observed_mz <- base::paste(
                                        temp_mz[match(sample_list,
                                                      temp_mz$Raw.file), 2],
                                        collapse=";"
                                    )
                                    temp_score <- stats::aggregate(sub$Score,
                                                                   by=lst,
                                                                   FUN=mean,
                                                                   na.rm=TRUE)
                                    # Order all observed scores according to
                                    # sample_list
                                    Observed_score <- base::paste(
                                        temp_score[match(sample_list,
                                                         temp_score$Raw.file),
                                                   2],
                                        collapse=";"
                                    )
                                    lst <- as.list(c(name, sequence, protein,
                                                    msmsscan,
                                                    base::paste(rownames(sub),
                                                                collapse=","),
                                                    base::paste(scores,
                                                                collapse=";"),
                                                    base::paste(num_matches,
                                                                collapse=";"),
                                                    Modifications,
                                                    Observed_RT,
                                                    Observed_mz,
                                                    Observed_score))
                                    data.table::set(
                                        x=features,
                                        i=as.integer(count_features),
                                        j=as.integer(c(1, 4, 5, 6, 11, 15, 16,
                                                       17, 19, 20, 21)),
                                        value=lst
                                    )

                                    mzb1 <- median_m.z + borders_m.z[1]
                                    mzb2 <- median_m.z + borders_m.z[2]
                                    rtb1 <- median_RT + borders_RT[1]
                                    rtb2 <- median_RT + borders_RT[2]
                                    lst <- as.list(c(as.numeric(median_m.z),
                                                     as.numeric(median_RT),
                                                     as.numeric(mzb1),
                                                     as.numeric(mzb2),
                                                     as.numeric(rtb1),
                                                     as.numeric(rtb2),
                                                     as.numeric(nrow(sub)),
                                                     as.numeric(median_mass),
                                                     as.numeric(Charge),
                                                     RT_length))
                                    data.table::set(
                                        x=features,
                                        i=as.integer(count_features),
                                        j=as.integer(c(2, 3, 7, 8, 9, 10, 12,
                                                       13, 14, 18)),
                                        value=lst)



                                    # Now mark matched features
                                    data.table::set(temp_data,
                                                    i=as.integer(rownames(sub)),
                                                    j=1L, value=1)
                                    # If RT_window had to be increased, reset
                                    # borders_RT
                                    if(too_large_RT_variation == TRUE){
                                        too_large_RT_variation <- FALSE
                                        borders_RT <- borders_RT_save
                                    }
                                }
                            }
                        }
                    }
                }

                updatecounter <- updatecounter + 1

                if(updatecounter >= 50){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))

                    updatecounter <- 0
                    label <- base::paste(round(i / max * 100, 0), " % done (",
                                         i, "/", max, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }
        }

        if(MassSpec_mode == "TIMSToF"){
            for(i in 1:length(unique_peptides)){
                lidx <- (last_index - 1 + num_rows_x_pept_seq[i, 2])
                ord <- order_by_Sequence[last_index:lidx]
                sub_peptide <- allpeptides[ord, ]
                last_index <- last_index + num_rows_x_pept_seq[i, 2]
                loc <- which(sub_peptide$Sequence == unique_peptides[i])
                sub_peptide <- sub_peptide[loc, ]

                for(c in unique(sub_peptide$Charge)){
                    for(m in unique(sub_peptide$Modifications)){
                        c1 <- sub_peptide$Charge == c
                        c2 <- sub_peptide$Modifications == m
                        sub <- sub_peptide[which(c1 & c2), ]

                        if(nrow(sub) > 0){
                            sub <- sub[order(sub$Intensity, decreasing=TRUE), ]

                            # Only use features close to median feature over all
                            # samples for this peptide
                            median_m.z <- stats::median(sub$m.z, na.rm=TRUE)
                            median_RT <- stats::median(sub$Retention.time,
                                                       na.rm=TRUE)
                            median_IM <- stats::median(sub$inverse_K0,
                                                       na.rm=TRUE)

                            c1 <- (abs(sub$inverse_K0 - median_IM)
                                   < borders_IM[2])
                            c3 <- abs(sub$m.z - median_m.z) < borders_m.z[2]
                            c2 <- (abs(sub$Retention.time - median_RT)
                                   < borders_RT[2])
                            selection <- which(c1 & c2 & c3)
                            sub <- sub[selection, ]

                            if(nrow(sub) > 0){
                                # Check if peptide was sequenced several times.
                                # If this is the case, use feature closest to
                                # median
                                if(any(duplicated(sub$Raw.file))){
                                    sub <- sub[!duplicated(sub$Raw.file), ]
                                }

                                median_m.z <- stats::median(sub$m.z, na.rm=TRUE)
                                median_RT <- stats::median(sub$Retention.time,
                                                           na.rm=TRUE)
                                median_IM <- stats::median(sub$inverse_K0,
                                                           na.rm=TRUE)
                                rt <- sub$Retention.time
                                RT_max_min <- (max(rt, na.rm=TRUE)
                                               - min(rt, na.rm=TRUE))
                                too_large_RT_variation <- FALSE

                                # Get subset of allpeptides based on charge and
                                # mz window indexing
                                x <- allpeptides_frag_indices_per_mz_window[[c]]
                                start_ind <- which(x[, 1] > median_m.z
                                                       + borders_m.z[1])[1] - 1
                                end_ind <- which(x[, 2] > median_m.z
                                                     + borders_m.z[2])[1] + 1

                                if(is.na(end_ind) | end_ind > nrow(x)){
                                    end_ind <- nrow(x)
                                }

                                c1 <- is.na(x[start_ind, 3])
                                c2 <- !is.na(x[start_ind + 1, 3])

                                # XXX: Double check this since indexes are being
                                # reassigned on each conditional
                                if(c1 & c2){
                                    start_ind <- x[start_ind + 1, 3]
                                    end_ind_temp <- x[end_ind, 4]

                                    if(is.na(end_ind_temp)){
                                        aux <- x[(end_ind + 1):nrow(x), 3]
                                        end_ind <- na.omit(aux)[1]
                                    }

                                    else{
                                        end_ind <- end_ind_temp
                                    }
                                }

                                c1 <- !is.na(x[start_ind, 3])
                                c2 <- is.na(x[end_ind, 4])

                                else if(c1 & c2){
                                    end_ind <- x[start_ind, 4]
                                    start_ind <- x[start_ind, 3]
                                }

                                c1 <- is.na(x[start_ind, 3])
                                c2 <- !is.na(x[end_ind, 3])

                                else if(c1 & c2){
                                    start_ind <- x[end_ind, 3]
                                    end_ind <- x[end_ind, 4]
                                }

                                else{
                                    start_ind <- x[start_ind, 3]
                                    end_ind <- x[end_ind, 4]
                                }

                                loc <- start_ind:end_ind
                                temp_sub <- allpeptides_frag[[c]][loc, ]
                                tmz <- temp_sub$m.z
                                c1 <- tmz >= median_m.z + borders_m.z[1]
                                c2 <- tmz <= median_m.z + borders_m.z[2]
                                trt <- temp_sub$Retention.time
                                c3 <- trt >= median_RT + borders_RT[1]
                                c4 <- trt <= median_RT + borders_RT[2]
                                tik0 <- temp_sub$inverse_K0
                                c5 <- tik0 >= median_IM + borders_IM[1]
                                c6 <- tik0 <= median_IM + borders_IM[2]
                                c7 <- temp_sub$Sequence == unique_peptides[i]
                                c8 <- temp_sub$Charge == c
                                c9 <- temp_sub$Modifications == m
                                sel <- which(c1 & c2 & c3 & c4 & c5 & c6 |
                                             c7 & c8 & c9)
                                sub <- temp_sub[sel, ]

                                sub <- sub[order(sub$Intensity, sub$Sequence,
                                                 decreasing=TRUE), ]

                                # Check if selected MaxQ features in expected
                                # RT, IM and mz window contains outliers (too
                                # high deviation in RT, IM or mz compared to all
                                # other features)
                                box_stats_rt <- grDevices::boxplot.stats(
                                    sub$Retention.time
                                )$stats
                                box_stats_mz <- grDevices::boxplot.stats(
                                    sub$m.z
                                )$stats
                                box_stats_im <- grDevices::boxplot.stats(
                                    sub$inverse_K0
                                )$stats

                                c1 <- sub$Retention.time < box_stats_rt[2]
                                c2 <- sub$Retention.time > box_stats_rt[4]
                                c3 <- sub$m.z < box_stats_mz[2]
                                c4 <- sub$m.z > box_stats_mz[4]
                                c5 <- sub$inverse_K0 < box_stats_im[2]
                                c6 <- sub$inverse_K0 > box_stats_im[4]
                                outlier <- which(c1 | c2 | c3 | c4 | c5 | c6)
                                # Features with sequence identification will be
                                # never regarded as outlier
                                c <- sub$Sequence[outlier] == unique_peptides[i]

                                if(any(c)){
                                    outlier <- outlier[-which(c)]
                                }

                                # Any outlier detected? exclude them
                                if(length(outlier) > 0){
                                    sub <- sub[-outlier, ]
                                }

                                # Check if the fraction of features with not
                                # currently searched sequence is > 25 %
                                # If yes, keep all sequences and indicate this
                                # with multiple sequences for this feature
                                # otherwise remove other sequences
                                loc <- which(sub$Sequence != " "
                                             & sub$Sequence != "")
                                count_seqs <- plyr::count(sub$Sequence[loc])
                                loc <- which(count_seqs$x == unique_peptides[i])
                                ff <- 0.8 * sum(count_seqs$freq)

                                if(count_seqs$freq[loc] >= ff){
                                    c1 <- sub$Sequence == unique_peptides[i]
                                    c2 <- sub$Sequence == " "
                                    c3 <- sub$Sequence == ""
                                    sub <- sub[which(c1 | c2 | c3), ]
                                }

                                if(any(duplicated(sub$Raw.file))){
                                    sub <- sub[!duplicated(sub$Raw.file), ]
                                }
                                # Variation in RT for same peptide sequence
                                # between samples is too high
                                if(nrow(sub) == 0){
                                    borders_RT_save <- borders_RT
                                    too_large_RT_variation <- TRUE
                                    borders_RT <- c(-RT_max_min / 2,
                                                    RT_max_min / 2)
                                    tmz <- temp_sub$m.z
                                    c1 <- tmz >= median_m.z + borders_m.z[1]
                                    c2 <- tmz <= median_m.z + borders_m.z[2]
                                    trt <- temp_sub$Retention.time
                                    c3 <- trt >= median_RT + borders_RT[1]
                                    c4 <- trt <= median_RT + borders_RT[2]
                                    tik0 <- temp_sub$inverse_K0
                                    c5 <- tik0 >= median_IM + borders_IM[1]
                                    c6 <- tik0 <= median_IM + borders_IM[2]
                                    loc <- which(c1 & c2 & c3 & c4 & c5 & c6)
                                    sub <- temp_sub[loc, ]

                                    # Still no match
                                    if(nrow(sub) == 0){
                                        too_large_RT_variation <- FALSE
                                        borders_RT <- borders_RT_save
                                    }
                                }

                                if(nrow(sub) > 0){
                                    count_features <- count_features + 1

                                    # New feature - give a new name
                                    name <- base::paste("Feature",
                                                        count_features, sep="_")
                                    nchr <- nchar(as.character(sub$Sequence))
                                    sequences_relevant <- which(nchr > 1)

                                    if(length(sequences_relevant) > 0){
                                        lst <- c(sub_peptide$Sequence[1],
                                                 sub[sequences_relevant,
                                                     "Sequence"])
                                        sequence <- base::paste(unique(lst),
                                                                collapse=";")
                                    }

                                    else{
                                        sequence <- sub_peptide$Sequence[1]
                                    }

                                    nchr <- nchar(as.character(sub$Proteins))
                                    proteins_relevant <- which(nchr > 1)

                                    if(length(proteins_relevant) > 0){
                                        lst <- c(sub_peptide$Proteins[1],
                                                 sub[proteins_relevant,
                                                     "Proteins"])
                                        protein <- base::paste(unique(lst),
                                                               collapse=";")
                                    }

                                    else{
                                        protein <- sub_peptide$Proteins[1]
                                    }

                                    sn <- sub$MSMS.Scan.Numbers
                                    nchr <- nchar(as.character(sn))
                                    aux <- which(nchr > 1)

                                    if(length(aux) > 0){
                                        tmp <- base::paste(
                                            sub[aux, "Raw.file"],
                                            "_msms_",
                                            sub[aux, "MSMS.Scan.Numbers"],
                                            sep=""
                                        )
                                        msmsscan <- base::paste(unique(tmp),
                                                                collapse=";")
                                    }

                                    else{
                                        msmsscan <- ""
                                    }

                                    median_mass <- stats::median(sub$Mass,
                                                                 na.rm=TRUE)
                                    Charge <- c

                                    scores <- NULL
                                    num_matches <- NULL
                                    seqs <- stringr::str_split(sequence, ";")

                                    for(s in unique(unlist(seqs))){
                                        loc <- which(sub$Sequence == s)
                                        scores <- append(scores,
                                                         mean(sub$Score[loc],
                                                              na.rm=TRUE))
                                        num_matches <- append(num_matches,
                                                              nrow(sub[loc, ]))
                                    }

                                    Modifications <- ifelse(m != "Unmodified",
                                                            m, "")
                                    RT_length <- max(sub$Retention.Length,
                                                     na.rm=TRUE)
                                    lst <- list(Raw.file=sub$Raw.file)
                                    temp_RT <- stats::aggregate(
                                        sub$Retention.time,
                                        by=lst,
                                        FUN=mean,
                                        na.rm=TRUE
                                    )
                                    # Order all observed RT according to
                                    # sample_list
                                    Observed_RT <- base::paste(
                                        temp_RT[match(sample_list,
                                                      temp_RT$Raw.file), 2],
                                        collapse=";"
                                    )
                                    temp_mz <- stats::aggregate(sub$m.z, by=lst,
                                                                FUN=mean,
                                                                na.rm=TRUE)
                                    # Order all observed mz according to
                                    # sample_list
                                    Observed_mz <- base::paste(
                                        temp_mz[match(sample_list,
                                                      temp_mz$Raw.file), 2],
                                        collapse=";"
                                    )
                                    temp_score <- stats::aggregate(sub$Score,
                                                                   by=lst,
                                                                   FUN=mean,
                                                                   na.rm=TRUE)
                                    # Order all observed scores according to
                                    # sample_list
                                    Observed_score <- base::paste(
                                        temp_score[match(sample_list,
                                                         temp_score$Raw.file),
                                                   2],
                                        collapse=";"
                                    )
                                    IM_length <- max(sub$inverse_K0_length,
                                                     na.rm=TRUE)
                                    temp_IM <- stats::aggregate(sub$inverse_K0,
                                                                by=lst,
                                                                FUN=mean,
                                                                na.rm=TRUE)
                                    # Order all observed IM according to
                                    # sample_list
                                    Observed_IM <- base::paste(
                                        temp_IM[match(sample_list,
                                                      temp_IM$Raw.file), 2],
                                        collapse=";"
                                    )
                                    lst <- as.list(c(name, sequence, protein,
                                                     msmsscan,
                                                     base::paste(rownames(sub),
                                                                 collapse=","),
                                                    base::paste(scores,
                                                                collapse=";"),
                                                    base::paste(num_matches,
                                                                collapse=";"),
                                                    Modifications, Observed_RT,
                                                    Observed_mz, Observed_score,
                                                    Observed_IM))
                                    data.table::set(
                                        x=features,
                                        i=as.integer(count_features),
                                        j=as.integer(c(1, 4, 5, 6, 11, 15, 16,
                                                       17, 19, 20, 21, 26)),
                                        value=lst
                                    )
                                    lst <- as.list(c(
                                        as.numeric(median_m.z),
                                        as.numeric(median_RT),
                                        as.numeric(median_m.z + borders_m.z[1]),
                                        as.numeric(median_m.z + borders_m.z[2]),
                                        as.numeric(median_RT + borders_RT[1]),
                                        as.numeric(median_RT + borders_RT[2]),
                                        as.numeric(nrow(sub)),
                                        as.numeric(median_mass),
                                        as.numeric(Charge),
                                        RT_length,
                                        median_IM,
                                        IM_length,
                                        as.numeric(median_IM + borders_IM[1]),
                                        as.numeric(median_IM + borders_IM[2])
                                    ))
                                    data.table::set(
                                        x=features,
                                        i=as.integer(count_features),
                                        j=as.integer(c(2, 3, 7, 8, 9, 10, 12,
                                                       13, 14, 18, 22, 23, 24,
                                                       25)),
                                        value = lst
                                    )

                                    # Now mark matched features
                                    data.table::set(temp_data,
                                                    i=as.integer(rownames(sub)),
                                                    j=1L, value=1)

                                    # If RT_window had to be increased, reset
                                    # borders_RT
                                    if(too_large_RT_variation == TRUE){
                                        too_large_RT_variation <- FALSE
                                        borders_RT <- borders_RT_save
                                    }
                                }
                            }
                        }
                    }
                }

                updatecounter <- updatecounter + 1

                if(updatecounter >= 50){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))
                    updatecounter <- 0
                    label <- base::paste(round(i / max * 100, 0), " % done (",
                                         i, "/", max, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }
        }

        close(pb)

        features <- features[1:count_features, ]
        borders_RT <- borders_RT_use_save
        fname <- base::paste("Temporary_files/features_aligned_step_1",
                             output_file_names_add, ".RData", sep="")

        if(MassSpec_mode == "Orbitrap"){
            save(temp_data, features, allpeptides, borders_m.z, borders_RT,
                 windows, file=fname)
        }

        if(MassSpec_mode == "TIMSToF"){
            save(temp_data, features, allpeptides, borders_m.z, borders_RT,
                 borders_IM, windows, file=fname)
        }

        # Prepare columns to indicate number of charges and with which other
        # features the respective feature was merged
        features$num_diff_charges <- NA
        features$num_diff_charges <- as.numeric(features$num_diff_charges)
        features$merged_with <- ""
        ncolumns <- ncol(features)

        # Collapse features with same m/z and RT parameters and count number of
        # different charge states for same feature
        borders_mass <- c(-feature_mass_deviation_collapse,
                          feature_mass_deviation_collapse)

        max <- nrow(features)
        pb <- tcltk::tkProgressBar(title="Merge features",
                                   label=base::paste(round(0 / max * 100, 0),
                                                     "% done"), min=0, max=max,
                                   width=300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0

        if(MassSpec_mode == "Orbitrap"){
            for(i in 1:nrow(features)){
                if(is.na(features[i,"num_diff_charges"])){
                    c1 <- features$Mass >= features$Mass[i] + borders_mass[1]
                    c2 <- features$Mass <= features$Mass[i] + borders_mass[2]
                    c3 <- features$RT >= features$RT[i] + borders_RT[1]
                    c4 <- features$RT <= features$RT[i] + borders_RT[2]
                    selection <- which(c1 & c2 & c3 & c4)

                    if(length(selection) > 1){
                        # Merge same charges
                        for(c in unique(features[selection, ]$Charge)){
                            # More than one feature with same charge -> merge
                            loc <- which(features[selection, ]$Charge == c)

                            if(length(loc) > 1){
                                fsel <- features[selection[loc], ]
                                m.z_range_min <- min(fsel$m.z_range_min,
                                                     na.rm=TRUE)
                                m.z_range_max <- max(fsel$m.z_range_max,
                                                     na.rm=TRUE)
                                RT_range_min <- min(fsel$RT_range_min,
                                                    na.rm=TRUE)
                                RT_range_max <- max(fsel$RT_range_max,
                                                    na.rm=TRUE)
                                aux <- unlist(fsel$ion_indices)
                                ion_indices <- base::paste(aux, collapse=",")
                                count_ion_indices <- sum(fsel$count_ion_indices,
                                                         na.rm=TRUE)
                                mass <- mean(unlist(fsel$Mass), na.rm=TRUE)
                                m.z <- mean(c(m.z_range_min, m.z_range_max))
                                RT <- mean(c(RT_range_min, RT_range_max))
                                seqs <- stringr::str_split(fsel$Sequence, ";")
                                Sequence <- unique(unlist(seqs))

                                scores <- NULL
                                num_matches <- NULL
                                Modifications <- NULL

                                for(s in unique(Sequence)){
                                    c1 <- features[selection, "Charge"] == c
                                    c2 <- grepl(s, features[selection,
                                                            "Sequence"])
                                    fts <- features[selection[which(c1 & c2)], ]
                                    sc <- stringr::str_split(fts$mean_Scores,
                                                             ";")
                                    scores_temp <- as.numeric(unlist(sc))
                                    scores <- append(scores, mean(scores_temp,
                                                                  na.rm=TRUE))
                                    nm <- stringr::str_split(fts$num_matches,
                                                             ";")
                                    num_matches_temp <- as.numeric(unlist(nm))
                                    num_matches <- append(num_matches,
                                                          sum(num_matches_temp,
                                                              na.rm=TRUE))
                                    md <- stringr::str_split(fts$Modifications,
                                                             ";")
                                    Modifications_temp <- unlist(md)
                                    mds <- base::paste(Modifications_temp,
                                                       collapse=",")
                                    Modifications <- append(Modifications, mds)
                                }

                                scores <- base::paste(scores, collapse=",")
                                num_matches <- base::paste(num_matches,
                                                           collapse=",")
                                Modifications <- base::paste(Modifications,
                                                             collapse=",")

                                if(any(Sequence == "")){
                                    Sequence <- Sequence[-which(Sequence == "")]
                                }

                                Sequence <- base::paste(Sequence, collapse=";")
                                prt <- stringr::str_split(fsel$Protein, ";")
                                Protein <- unique(unlist(prt))

                                if(any(Protein == "")){
                                    Protein <- Protein[-which(Protein == "")]
                                }

                                Protein <- base::paste(Protein, collapse=";")
                                sn <- stringr::str_split(fsel$MSMS.Scan.Numbers,
                                                         ";")
                                MSMS.Scan.Numbers <- unique(unlist(sn))

                                if(any(MSMS.Scan.Numbers == "")){
                                    loc <- which(MSMS.Scan.Numbers == "")
                                    MSMS.Scan.Numbers <- MSMS.Scan.Numbers[-loc]
                                }

                                MSMS.Scan.Numbers <- base::paste(
                                    MSMS.Scan.Numbers,
                                    collapse=";"
                                )
                                Observed_RTs <- stringr::str_split(
                                    fsel$Observed_RT,
                                    pattern=";",
                                    simplify=TRUE
                                )
                                Observed_RTs <- apply(Observed_RTs, 2,
                                                      as.numeric)
                                Observed_RTs <- colMeans(Observed_RTs,
                                                         na.rm=TRUE)
                                Observed_RTs[is.na(Observed_RTs)] <- NA
                                Observed_RTs <- base::paste(Observed_RTs,
                                                            collapse=";")
                                Observed_mz <- stringr::str_split(
                                    fsel$Observed_mz,
                                    pattern=";",
                                    simplify=TRUE
                                )
                                Observed_mz <- apply(Observed_mz, 2, as.numeric)
                                Observed_mz <- colMeans(Observed_mz, na.rm=TRUE)
                                Observed_mz[is.na(Observed_mz)] <- NA
                                Observed_mz <- base::paste(Observed_mz,
                                                           collapse=";")
                                Observed_score <- stringr::str_split(
                                    fsel$Observed_score,
                                    pattern=";",
                                    simplify=TRUE
                                )
                                Observed_score <- apply(Observed_score, 2,
                                                        as.numeric)
                                Observed_score <- colMeans(Observed_score,
                                                           na.rm=TRUE)
                                Observed_score[is.na(Observed_score)] <- NA
                                Observed_score <- base::paste(Observed_score,
                                                              collapse=";")

                                js <- c(2, 3, 7, 8, 9, 10, 12, 13, 14)
                                data.table::set(
                                    x=features,
                                    i=selection[loc[1]],
                                    j=as.integer(js),
                                    value=as.list(c(m.z, RT, m.z_range_min,
                                                    m.z_range_max, RT_range_min,
                                                    RT_range_max,
                                                    count_ion_indices, mass, c))
                                )

                                js <- c(4, 5, 6, 11, 15, 16, 17, 19, 20, 21)
                                data.table::set(
                                    x=features,
                                    i=selection[loc[1]],
                                    j=as.integer(js),
                                    value=as.list(c(Sequence, Protein,
                                                    MSMS.Scan.Numbers,
                                                    ion_indices, scores,
                                                    num_matches, Modifications,
                                                    Observed_RTs, Observed_mz,
                                                    Observed_score))
                                )

                                val <- features$Feature_name[selection[loc[1]]]
                                data.table::set(
                                    x=features,
                                    i=selection[loc[-1]],
                                    j=as.integer(ncolumns),
                                    value=val
                                )
                            }
                        }

                        #Add number of different charges count
                        data.table::set(
                            x=features,
                            i=selection,
                            j=as.integer(ncolumns - 1),
                            value=length(unique(features[selection, ]$Charge))
                        )
                    }
                }

                updatecounter <- updatecounter + 1

                if(updatecounter >= 100){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))
                    updatecounter <- 0
                    label <- base::paste(round(i / max * 100, 0), " % done (",
                                         i, "/", max, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }
        }

        if(MassSpec_mode == "TIMSToF"){
            for(i in 1:nrow(features)){
                if(is.na(features[i, "num_diff_charges"])){
                    c1 <- features$Mass >= features$Mass[i] + borders_mass[1]
                    c2 <- features$Mass <= features$Mass[i] + borders_mass[2]
                    c3 <- features$RT >= features$RT[i] + borders_RT[1]
                    c4 <- features$RT <= features$RT[i] + borders_RT[2]
                    c5 <- features$Inv_K0 >= features$Inv_K0[i] + borders_IM[1]
                    c6 <- features$Inv_K0 <= features$Inv_K0[i] + borders_IM[2]
                    selection <- which(c1 & c2 & c3 & c4 & c5 & c6)

                    if(length(selection) > 1){
                        # Merge same charges
                        for(c in unique(features[selection, ]$Charge)){
                            # More than one feature with same charge -> merge
                            aux <- which(features[selection, ]$Charge == c)

                            if(length(aux) > 1){
                                fts <- features[selection[aux], ]
                                m.z_range_min <- min(fts$m.z_range_min,
                                                     na.rm=TRUE)
                                m.z_range_max <- max(fts$m.z_range_max,
                                                     na.rm=TRUE)
                                RT_range_min <- min(fts$RT_range_min,
                                                    na.rm=TRUE)
                                RT_range_max <- max(fts$RT_range_max,
                                                    na.rm=TRUE)
                                IM_range_min <- min(fts$Inv_K0_range_min,
                                                    na.rm=TRUE)
                                IM_range_max <- max(fts$Inv_K0_range_max,
                                                    na.rm=TRUE)
                                inds <- unlist(fts$ion_indices)
                                ion_indices <- base::paste(inds, collapse=",")
                                count_ion_indices <- sum(fts$count_ion_indices,
                                                         na.rm=TRUE)
                                mass <- mean(unlist(fts$Mass), na.rm=TRUE)
                                m.z <- mean(c(m.z_range_min, m.z_range_max))
                                RT <- mean(c(RT_range_min, RT_range_max))
                                IM <- mean(c(IM_range_min, IM_range_max))
                                sqs <- stringr::str_split(fts$Sequence, ";")
                                Sequence <- unique(unlist(sqs))

                                scores <- NULL
                                num_matches <- NULL
                                Modifications <- NULL

                                for(s in unique(Sequence)){
                                    c1 <- features[selection, "Charge"] == c
                                    c2 <- grepl(s, features[selection,
                                                            "Sequence"])
                                    loc <- selection[which(c1 & c2)]
                                    feats <- features[loc, ]
                                    ms <- stringr::str_split(feats$mean_Scores,
                                                             ";")
                                    scores_temp <- as.numeric(unlist(ms))
                                    scores <- append(scores, mean(scores_temp,
                                                                  na.rm=TRUE))
                                    nm <- stringr::str_split(feats$num_matches,
                                                             ";")
                                    num_matches_temp <- as.numeric(unlist(nm))
                                    num_matches <- append(num_matches,
                                                          sum(num_matches_temp,
                                                              na.rm=TRUE))
                                    m <- stringr::str_split(feats$Modifications,
                                                            ";")
                                    Modifications_temp <- unlist(m)
                                    mt <- base::paste(Modifications_temp,
                                                      collapse=",")
                                    Modifications <- append(Modifications, mt)
                                }

                                scores <- base::paste(scores, collapse=",")
                                num_matches <- base::paste(num_matches,
                                                           collapse=",")
                                Modifications <- base::paste(Modifications,
                                                             collapse=",")

                                if(any(Sequence == "")){
                                    Sequence <- Sequence[-which(Sequence == "")]
                                }

                                Sequence <- base::paste(Sequence, collapse=";")
                                prt <- stringr::str_split(fts$Protein, ";")
                                Protein <- unique(unlist(prt))

                                if(any(Protein == "")){
                                    Protein <- Protein[-which(Protein == "")]
                                }

                                Protein <- base::paste(Protein, collapse=";")
                                sn <- stringr::str_split(fts$MSMS.Scan.Numbers,
                                                         ";")
                                MSMS.Scan.Numbers <- unique(unlist(sn))

                                if(any(MSMS.Scan.Numbers == "")){
                                    loc <- which(MSMS.Scan.Numbers == "")
                                    MSMS.Scan.Numbers <- MSMS.Scan.Numbers[-loc]
                                }

                                MSMS.Scan.Numbers <- base::paste(
                                    MSMS.Scan.Numbers,
                                    collapse=";"
                                )

                                Observed_RTs <- stringr::str_split(
                                    fts$Observed_RT,
                                    pattern=";",
                                    simplify=TRUE
                                )
                                Observed_RTs <- apply(Observed_RTs, 2,
                                                      as.numeric)
                                Observed_RTs <- colMeans(Observed_RTs,
                                                         na.rm=TRUE)
                                Observed_RTs[is.na(Observed_RTs)] <- NA
                                Observed_RTs <- base::paste(Observed_RTs,
                                                            collapse=";")
                                Observed_mz <- stringr::str_split(
                                    fts$Observed_mz,
                                    pattern=";",
                                    simplify=TRUE
                                )
                                Observed_mz <- apply(Observed_mz, 2, as.numeric)
                                Observed_mz <- colMeans(Observed_mz, na.rm=TRUE)
                                Observed_mz[is.na(Observed_mz)] <- NA
                                Observed_mz <- base::paste(Observed_mz,
                                                           collapse=";")
                                Observed_score <- stringr::str_split(
                                    fts$Observed_score,
                                    pattern=";",
                                    simplify=TRUE
                                )
                                Observed_score <- apply(Observed_score, 2,
                                                        as.numeric)
                                Observed_score <- colMeans(Observed_score,
                                                           na.rm=TRUE)
                                Observed_score[is.na(Observed_score)] <- NA
                                Observed_score <- base::paste(Observed_score,
                                                              collapse=";")
                                Observed_IMs <- stringr::str_split(
                                    fts$Observed_IM,
                                    pattern=";",
                                    simplify=TRUE
                                )
                                Observed_IMs <- apply(Observed_IMs, 2,
                                                      as.numeric)
                                Observed_IMs <- colMeans(Observed_IMs,
                                                         na.rm=TRUE)
                                Observed_IMs[is.na(Observed_IMs)] <- NA
                                Observed_IMs <- base::paste(Observed_IMs,
                                                            collapse=";")
                                js <- as.integer(c(2, 3, 7, 8, 9, 10, 12, 13,
                                                   14, 22, 24, 25))
                                lst <- as.list(c(m.z, RT, m.z_range_min,
                                                 m.z_range_max, RT_range_min,
                                                 RT_range_max,
                                                 count_ion_indices, mass, c, IM,
                                                 IM_range_min, IM_range_max))
                                data.table::set(x=features, i=selection[aux[1]],
                                                j=js, value=lst)
                                js <- as.integer(c(4, 5, 6, 11, 15, 16, 17, 19,
                                                   20, 21, 26))
                                lst <- as.list(c(Sequence, Protein,
                                                 MSMS.Scan.Numbers, ion_indices,
                                                 scores, num_matches,
                                                 Modifications, Observed_RTs,
                                                 Observed_mz, Observed_score,
                                                 Observed_IMs))
                                data.table::set(x=features, i=selection[aux[1]],
                                                j=js, value=lst)
                                loc <- which(features[selection, ]$Charge == c)
                                val <- features$Feature_name[selection[loc[1]]]
                                data.table::set(x=features,
                                                i=selection[aux[-1]],
                                                j=as.integer(ncolumns),
                                                value=val)
                            }
                        }
                        # Add number of different charges count
                        val <- length(unique(features[selection, ]$Charge))
                        data.table::set(x=features, i=selection,
                                        j=as.integer(ncolumns - 1),
                                        value=val)
                    }
                }

                updatecounter <- updatecounter + 1
                if(updatecounter >= 100){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))
                    updatecounter <- 0
                    label <- base::paste(round(i / max * 100, 0), " % done (",
                                         i, "/", max, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }
        }

        close(pb)

        features$num_diff_charges[is.na(features$num_diff_charges)] <- 1
        features <- features[which(features$merged_with == ""), ]
        features <- base::as.data.frame(features)[, -ncolumns]
        # Now add mz and RT calibration information per feature
        mat <- matrix(ncol=length(sample_list), nrow=nrow(features), NA)
        RT_calibration_vals <- base::as.data.frame(mat)
        colnames(RT_calibration_vals) <- base::paste("RT_calibration",
                                                     sample_list, sep=" ")

        for(c in 1:ncol(RT_calibration_vals)){
            RT_calibration_vals[, c] <- as.numeric(RT_calibration_vals[, c])
        }

        mat <- matrix(ncol=length(sample_list), nrow=nrow(features), NA)
        mz_calibration_vals <- base::as.data.frame(mat)
        colnames(mz_calibration_vals) <- base::paste("mz_calibration",
                                                     sample_list, sep=" ")
        for(c in 1:ncol(mz_calibration_vals)){
            mz_calibration_vals[, c] <- as.numeric(mz_calibration_vals[, c])
        }

        if(MassSpec_mode == "TIMSToF"){
            mat <- matrix(ncol=length(sample_list), nrow=nrow(features), NA)
            IM_calibration_vals <- base::as.data.frame(mat)
            colnames(IM_calibration_vals) <- base::paste("IM_calibration",
                                                         sample_list, sep=" ")

            for(c in 1:ncol(IM_calibration_vals)){
                IM_calibration_vals[, c] <- as.numeric(IM_calibration_vals[, c])
            }
        }

        # Any of both calibrations should be done?
        if(any(c(RT_calibration, mz_calibration) == TRUE)){
            # First: start by taking already available calibration information
            peptide_features <- which(features$Sequence != "")
            pep_seq <- stringr::str_split(features$Sequence, ";", simplify=TRUE)
            evidence <- evidence[which(evidence$Raw.file %in% sample_list), ]

            # Evidence_peptides
            evidence_peptides <- unique(evidence$Sequence)
            len <- length(evidence_peptides)
            evidence_peptides_indices <- vector(mode="list", length=len)
            names(evidence_peptides_indices) <- evidence_peptides

            max <- length(evidence_peptides)
            label <- base::paste(round(0 / max * 100, 0), "% done")
            pb <- tcltk::tkProgressBar(title="Indexing peptides", label=label,
                                       min=0, max=max, width=300)
            start_time <- Sys.time()
            updatecounter <- 0
            time_require <- 0
            start <- 1

            for(p in evidence_peptides){
                i <- which(evidence_peptides == p)
                evidence_peptides_indices[[p]] <- which(evidence$Sequence == p)
                updatecounter <- updatecounter + 1

                if(updatecounter >= 100){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))
                    updatecounter <- 0
                    label <- base::paste(round(i / max * 100, 0), " % done (",
                                         i, "/", max, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }

            close(pb)

            # Grap already available information for RT and mz calibration from
            # MaxQ output
            max <- length(peptide_features)
            title <- "Preparing calibrations per feature"
            label <- base::paste(round(0 / max * 100, 0), "% done")
            pb <- tcltk::tkProgressBar(title=title, label=label, min=0, max=max,
                                       width=300)
            start_time <- Sys.time()
            updatecounter <- 0
            time_require <- 0
            counter <- 0

            for(i in peptide_features){
                cur_mods <- stringr::str_split(features$Modifications[i], ";|,",
                                               simplify=TRUE)
                cur_mods[cur_mods == ""] <- "Unmodified"
                cur_mods <- unique(cur_mods)
                indices <- NULL

                for(p in pep_seq[i, which(pep_seq[i, ] != "")]){
                    indices <- append(indices, evidence_peptides_indices[[p]])
                }

                evidence_sub <- evidence[indices, ]
                c1 <- evidence_sub$Charge == features$Charge[i]
                c2 <- evidence_sub$Modifications %in% cur_mods
                indices <- which(c1 & c2)

                if(length(indices) > 0){
                    if(RT_calibration == TRUE){
                        rts <- stringr::str_split(features$Observed_RT[i],
                                                  pattern=";", simplify=TRUE)
                        calc_RT <- as.numeric(rts) - features$RT[i]
                        data.table::set(RT_calibration_vals, as.integer(i),
                                        as.integer(1:length(sample_list)),
                                        value=as.list(c(calc_RT)))
                    }

                    if(IM_calibration == TRUE){
                        ims <- stringr::str_split(features$Observed_IM[i],
                                                  pattern=";", simplify=TRUE)
                        calc_IM <- as.numeric(ims) - features$Inv_K0[i]
                        data.table::set(IM_calibration_vals, as.integer(i),
                                        as.integer(1:length(sample_list)),
                                        value=as.list(c(calc_IM)))
                    }

                    if(mz_calibration == TRUE){
                        # Determine mz correction based on observed mz peaks
                        if(use_mz_at_max_int_for_correction == TRUE){
                            mzs <- stringr::str_split(features$Observed_mz[i],
                                                      pattern=";",
                                                      simplify=TRUE)
                            calc_mz <- as.numeric(mzs) - features$m.z[i]
                            data.table::set(mz_calibration_vals, as.integer(i),
                                            as.integer(1:length(sample_list)),
                                            value=as.list(c(calc_mz)))
                        }
                        # Use MaxQ correction
                        else{
                            col <- "Uncalibrated...Calibrated.m.z..Da."
                            calc_mz <- stats::aggregate(
                                evidence_sub[indices, col],
                                by=list(Sample=evidence_sub$Raw.file[indices]),
                                FUN="mean",
                                na.rm=TRUE
                            )
                            data.table::set(mz_calibration_vals, as.integer(i),
                                            as.integer(match(calc_mz$Sample,
                                                             sample_list)),
                                            value=as.list(c(calc_mz$x)))
                        }
                    }
                }

                counter <- counter + 1
                updatecounter <- updatecounter + 1

                if(updatecounter >= 100){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- ((time_elapsed / (counter / max))
                                     * (1 - (counter / max)))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))
                    updatecounter <- 0
                    label <- base::paste(round(counter / max * 100, 0),
                                         " % done (", counter, "/", max,
                                         ", Time require: ", time_require, ")",
                                         sep="")
                    tcltk::setTkProgressBar(pb, counter, label=label)
                }
            }

            close(pb)
            RT_calibration_vals_save <- RT_calibration_vals
            c1 <- mz_calibration == TRUE
            c2 <- use_mz_at_max_int_for_correction == FALSE

            if(c1 & c2){
                # Train random forest model for predicting mz calibration
                # Train a model per sample
                models <- list()
                trainsets <- list()
                evalsets <- list()
                ucmd <- "Uncalibrated...Calibrated.m.z..Da."

                if(any(colnames(evidence) == "Resolution")){
                    usecols <- c("Raw.file", "Retention.time", "m.z", "Charge",
                                 "Resolution", ucmd)
                    ev <- evidence[, usecols]
                    temp_data <- evidence[which(rowSums(is.na(ev)) == 0),
                                          usecols]
                    temp_data$Resolution <- 0
                }

                else{
                    usecols <- c("Raw.file", "Retention.time", "m.z", "Charge",
                                 ucmd)
                    ev <- evidence[, usecols]
                    temp_data <- evidence[which(rowSums(is.na(ev)) == 0),
                                          usecols]
                    temp_data$Resolution <- 0
                }

                print(paste0(Sys.time(), " Train RF models for mz-corrections"))

                for(s in sample_list){
                    selec_rows <- which(temp_data$Raw.file == s)
                    temp_data2 <- temp_data[selec_rows, ]
                    random_selection <- sample(1:nrow(temp_data2),
                                               size=0.8 * nrow(temp_data2),
                                               replace=FALSE)
                    trainset <- temp_data2[random_selection, ]
                    evalset <- temp_data2[-random_selection, ]
                    unq <- unique(trainset[, ucmd])

                    if(length(unq) > 5){
                        tsetx <- trainset[, c("Retention.time", "m.z", "Charge",
                                              "Resolution")]
                        tsety <- trainset[, ucmd]
                        model <- randomForest::randomForest(tsetx, tsety,
                                                            importance=TRUE,
                                                            ntree=500, mtry=4,
                                                            do.trace=FALSE,
                                                            nodesize=100)
                        evalset$predicted <- stats::predict(model, evalset,
                                                            type="response")
                        trainset$predicted <- stats::predict(model, trainset,
                                                             type="response")
                        models[[s]] <- model
                        trainsets[[s]] <- trainset
                        evalsets[[s]] <- evalset

                        title <- base::paste(s, "- Train set (80% of data)")
                        graphics::plot(
                            trainset$predicted,
                            trainset$Uncalibrated...Calibrated.m.z..Da.,
                            main=title,
                            xlab="Predicted m/z calibration",
                            ylab="MaxQ determined m/z calibration"
                        )
                        aux <- trainset$Uncalibrated...Calibrated.m.z..Da.
                        fit <- stats::lm(aux~trainset$predicted)
                        graphics::abline(fit)
                        posx <- (graphics::par("usr")[1]
                                 + (graphics::par("usr")[2]
                                    - graphics::par("usr")[1]) * 0.8)
                        posy <- (graphics::par("usr")[3]
                                 + (graphics::par("usr")[4]
                                    - graphics::par("usr")[3]) * 0.2)
                        r2t <- round(as.numeric(summary(fit)[8]), digits=2)
                        graphics::text(posx, posy, base::paste("R^2:", r2t))
                        fit_train <- fit

                        title <- base::paste(s,
                                             "- Validation set (20% of data)")
                        graphics::plot(
                            evalset$predicted,
                            evalset$Uncalibrated...Calibrated.m.z..Da.,
                            main=title,
                            xlab="Predicted m/z calibration",
                            ylab="MaxQ determined m/z calibration"
                        )
                        aux <- evalset$Uncalibrated...Calibrated.m.z..Da.
                        fit <- stats::lm(aux~evalset$predicted)
                        graphics::abline(fit)
                        posx <- (graphics::par("usr")[1]
                                 + (graphics::par("usr")[2]
                                    - graphics::par("usr")[1]) * 0.8)
                        posy <- (graphics::par("usr")[3]
                                 + (graphics::par("usr")[4]
                                    - graphics::par("usr")[3]) * 0.2)
                        r2e <- round(as.numeric(summary(fit)[8]), digits=2)
                        graphics::text(posx,posy,base::paste("R^2:", r2e))
                        fit_eval <- fit

                        print(base::paste(s, ": R^2 train-set=", r2t,
                                          ", R^2 eval-set=", r2e, sep=""))
                    }
                }

                QC_data[["mz_calibration_models"]] <- models
                QC_data[["mz_calibration_train_sets"]] <- trainsets
                QC_data[["mz_calibration_evaluation_sets"]] <- evalsets

                # Now predict mz for missing mzcalibrations of features

                # Get ion indices per feature
                mat <- matrix(ncol=6, nrow=nrow(allpeptides))
                all_ion_indices <- base::as.data.frame(mat)
                colnames(all_ion_indices) <- c("Feature", "Raw.file",
                                               "Retention.time", "m.z",
                                               "Charge", "Resolution")
                all_ion_indices$Feature <- as.numeric(all_ion_indices$Feature)
                all_ion_indices$Raw.file <- as.character(
                    all_ion_indices$Raw.file
                )
                all_ion_indices$Retention.time <- as.numeric(
                    all_ion_indices$Retention.time
                )
                all_ion_indices$m.z <- as.numeric(all_ion_indices$m.z)
                all_ion_indices$Charge <- as.integer(all_ion_indices$Charge)
                all_ion_indices$Resolution <- as.numeric(
                    all_ion_indices$Resolution
                )
                cur_index <- 1

                max <- nrow(features)
                title <- "Prepare features for m/z-calibration prediction by RF"
                pb <- tcltk::tkProgressBar(
                    title=title,
                    label=base::paste(round(0 / max * 100, 0), "% done"),
                    min=0,
                    max=max,
                    width=300
                )
                start_time <- Sys.time()
                updatecounter <- 0
                time_require <- 0

                for(i in 1:nrow(features)){
                    # Start with extracting ion_indices per feature
                    ii <- stringr::str_split(features$ion_indices[i], ",",
                                             simplify=TRUE)
                    cur_ion_indices <- as.numeric(ii)
                    usei <- as.integer(cur_index:(cur_index - 1
                                                  + length(cur_ion_indices)))
                    data.table::set(all_ion_indices, i=usei, 1L, i)
                    data.table::set(all_ion_indices, i=usei, 2L,
                                    allpeptides[cur_ion_indices, "Raw.file"])
                    val <- allpeptides[cur_ion_indices, c("Retention.time",
                                                          "m.z", "Charge",
                                                          "Resolution")]
                    data.table::set(all_ion_indices, i=usei, as.integer(3:6),
                                    value=as.list(val))

                    cur_index <- cur_index + length(cur_ion_indices)

                    updatecounter <- updatecounter + 1

                    if(updatecounter >= 100){
                        time_elapsed <- difftime(Sys.time(), start_time,
                                                 units="secs")
                        time_require <- ((time_elapsed / (i / max))
                                         * (1 - (i / max)))
                        td <- lubridate::seconds_to_period(time_require)
                        time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                                lubridate::minute(td),
                                                round(lubridate::second(td),
                                                      digits=0))
                        updatecounter <- 0
                        label <- base::paste(round(i / max * 100, 0),
                                             " % done (", i, "/", max,
                                             ", Time require: ", time_require,
                                             ")", sep="")
                        tcltk::setTkProgressBar(pb, i, label=label)
                    }
                }

                close(pb)
                all_ion_indices <- all_ion_indices[1:(cur_index - 1), ]
                all_ion_indices$Resolution <- 0
                raw_files <- sample_list

                # Get mean ion propperties per feature to predict m/z
                # calibration per feature for each sample
                median_feature_properties <- all_ion_indices[, -2]

                # Determine mean values per feature
                mfp <- median_feature_properties
                lst <- list(Feature=mfp$Feature)
                median_feature_properties <- stats::aggregate(mfp[, -1], by=lst,
                                                              FUN="median",
                                                              na.rm=TRUE)
                QC_data[["mz_calibration_median_feature_properties"]] <- mfp
                max <- ncol(mz_calibration_vals)
                title <- "Predict mz calibrations using RandomForest models"
                label <- base::paste(round(0 / max * 100, 0), "% done")
                pb <- tcltk::tkProgressBar(title=title, label=label, min=0,
                                           max=max, width=300)
                start_time <- Sys.time()
                updatecounter <- 0
                time_require <- 0

                for(c in 1:ncol(mz_calibration_vals)){
                    select_model <- which(names(models) == raw_files[c])

                    if(length(select_model) > 0){
                        ina <- is.na(mz_calibration_vals[, c])
                        features_select <- which(ina)
                        loc <- which(median_feature_properties$Feature
                                     %in% features_select)
                        temp_data <- median_feature_properties[loc, ]
                        rs <- rowSums(is.na(temp_data[, c("Retention.time",
                                                          "m.z", "Charge",
                                                          "Resolution")]))
                        temp_data <- temp_data[which(rs == 0),]

                        if(nrow(temp_data) > 0){
                            prediction <- stats::predict(models[[select_model]],
                                                         temp_data[, -1],
                                                         type="response")
                            data.table::set(mz_calibration_vals,
                                            as.integer(temp_data$Feature),
                                            as.integer(c), value=prediction)
                        }
                    }

                    updatecounter <- updatecounter + 1

                    if(updatecounter >= 1){
                        time_elapsed <- difftime(Sys.time(), start_time,
                                                 units="secs")
                        time_require <- ((time_elapsed / (c / max))
                                         * (1 - (c / max)))
                        td <- lubridate::seconds_to_period(time_require)
                        time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                                lubridate::minute(td),
                                                round(lubridate::second(td),
                                                      digits=0))

                        updatecounter <- 0
                        label <- base::paste(round(c / max * 100, 0),
                                             " % done (", c, "/", max,
                                             ", Time require: ", time_require,
                                             ")", sep="")
                        tcltk::setTkProgressBar(pb, c, label=label)
                    }
                }

                close(pb)
            }

            # If there are still missing values just fill with 0 or expected
            # SILAC shifts
            if(any(is.na(mz_calibration_vals))){
                mz_calibration_vals[is.na(mz_calibration_vals)] <- 0
            }

            if(RT_calibration == TRUE){
                graphics::plot(1, type="n", axes=FALSE, xlab="", ylab="")
                graphics::text(1, 1,
                               "RT GAM fitting - First round of estimation")

                RT_alignment_GAM_models <- list()
                # Fit GAM and predict per sample.
                for(c in 1:ncol(RT_calibration_vals)){
                    x <- features$RT[which(RT_calibration_vals[, c] != 0)]
                    loc <- which(RT_calibration_vals[, c] != 0)
                    y <- RT_calibration_vals[, c][loc]
                    # Require at least 500 data points to perform fitting
                    if(length(x) > 10){
                        # Try to fit an average generalised additive model to
                        # determine a RT dependent calibration curve
                        gam <- mgcv::gam(y ~ s(x), method="REML")
                        RT_alignment_GAM_models[[c]] <- gam
                        x_pred <- seq(min(features$RT, na.rm=TRUE),
                                      max(features$RT, na.rm=TRUE),
                                      length.out=nrow(features))
                        y_pred <- stats::predict(gam,
                                                 base::data.frame(x=x_pred))
                        lim_stats <- grDevices::boxplot.stats(y_pred)
                        delta_y <- lim_stats$stats[4] - lim_stats$stats[2]
                        ylim <- c(-max(abs(c(lim_stats$stats[2],
                                             lim_stats$stats[4]))) - delta_y,
                                  max(abs(c(lim_stats$stats[2],
                                            lim_stats$stats[4]))) + delta_y)
                        if(ylim[2] < 2){
                            ylim <- c(-2, 2)
                        }

                        graphics::smoothScatter(x, y,
                                                ylab="Observed RT calibration",
                                                main=sample_list[c],
                                                xlab="RT [min]", ylim=ylim)
                        graphics::lines(x_pred, y_pred, col="red")
                        graphics::legend("topright", legend="GAM", lty=1,
                                         col="red")

                        # Predict RT correction for all features in current
                        # sample
                        df <- base::data.frame(x=features$RT)
                        y_pred <- stats::predict(gam, df)

                        # Use GAM to predict RT correction for missing features
                        # in current sample
                        selection <- which(is.na(RT_calibration_vals[, c]))
                        RT_calibration_vals[, c][selection] <- y_pred[selection]
                    }
                    # Not enough observations ... skip
                    else{
                        msg <- paste0("- Not enough peptide observations for ",
                                      "RT-GAM fitting...")
                        print(base::paste(sample_list[c], msg))
                    }
                }

                col <- "RT_calibration_GAM_models"
                QC_data[[col]] <- RT_alignment_GAM_models

                # Perform a second round of RT calibration
                # Currently it is possible that several features are observed
                # only in a few samples hence observed median RT is driven by
                # these few samples based on GAM models adjust median feature RT
                # and rerun fitting
                expected_raw_RTs <- features$RT + RT_calibration_vals
                expected_raw_RTs[is.na(expected_raw_RTs)] <- 0

                delta_rt_temp <- features$RT_range_max - features$RT
                mat <- as.matrix(expected_raw_RTs)
                features$RT <- matrixStats::rowMedians(mat, na.rm=TRUE)
                features$RT_range_min <- features$RT - delta_rt_temp
                features$RT_range_max <- features$RT + delta_rt_temp

                # Restart RT calibration
                # Grap known RT once again
                RT_calibration_vals <- RT_calibration_vals_save

                graphics::plot(1, type="n", axes=FALSE, xlab="", ylab="")
                graphics::text(1, 1,
                               "RT GAM fitting - Second round of estimation")
                RT_alignment_GAM_models <- list()

                # Fit GAM and predict per sample.
                for(c in 1:ncol(RT_calibration_vals)){
                    x <- features$RT[which(RT_calibration_vals[, c] != 0)]
                    loc <- which(RT_calibration_vals[, c] != 0)
                    y <- RT_calibration_vals[, c][loc]
                    # Require at least 500 data points to perform fitting
                    if(length(x) > 10){
                        # Try to fit an average generalised additive model to
                        # determine a RT dependent calibration curve
                        gam <- mgcv::gam(y ~ s(x), method="REML")
                        RT_alignment_GAM_models[[c]] <- gam
                        x_pred <- seq(min(features$RT, na.rm=TRUE),
                                      max(features$RT, na.rm=TRUE),
                                      length.out=nrow(features))
                        df <- base::data.frame(x=x_pred)
                        y_pred <- stats::predict(gam, df)
                        lim_stats <- grDevices::boxplot.stats(y_pred)
                        delta_y <- lim_stats$stats[4] - lim_stats$stats[2]
                        ylim <- c(-max(abs(c(lim_stats$stats[2],
                                             lim_stats$stats[4]))) - delta_y,
                                  max(abs(c(lim_stats$stats[2],
                                            lim_stats$stats[4]))) + delta_y)

                        if(ylim[2] < 2){
                            ylim <- c(-2, 2)
                        }

                        graphics::smoothScatter(x, y,
                                                ylab="Observed RT calibration",
                                                main=sample_list[c],
                                                xlab="RT [min]", ylim=ylim)
                        graphics::lines(x_pred, y_pred, col="red")
                        graphics::legend("topright", legend="GAM", lty=1,
                                         col="red")

                        # Predict RT correction for all features in current
                        # sample
                        df <- base::data.frame(x=features$RT)
                        y_pred <- stats::predict(gam, df)

                        # Use GAM to predict RT correction for missing features
                        # in current sample
                        selection <- which(is.na(RT_calibration_vals[, c]))
                        RT_calibration_vals[, c][selection] <- y_pred[selection]
                    }

                    # Not enough observations ... skip
                    else{
                        msg <- paste0("- Not enough peptide observations for ",
                                      "RT-GAM fitting...")
                        print(base::paste(sample_list[c], msg))
                    }
                }

                col <- "RT_calibration_GAM_models"
                QC_data[[col]] <- RT_alignment_GAM_models
            }

            RT_calibration_vals[is.na(RT_calibration_vals)] <- 0

            if(MassSpec_mode == "TIMSToF" & IM_calibration == TRUE){
                # Use median of IM_calibration per sample for NAs
                for(c in 1:ncol(IM_calibration_vals)){
                    loc <- is.na(IM_calibration_vals[, c])
                    IM_calibration_vals[, c][loc] <- stats::median(
                        IM_calibration_vals[, c],
                        na.rm=TRUE
                    )
                }
            }
        }
        # No calibration should be done
        else{
            RT_calibration_vals[is.na(RT_calibration_vals)] <- 0
            mz_calibration_vals[is.na(mz_calibration_vals)] <- 0

            if(MassSpec_mode == "TIMSToF"){
                IM_calibration_vals[is.na(IM_calibration_vals)] <- 0
            }
        }

        features <- cbind(features, RT_calibration_vals, mz_calibration_vals)

        if(MassSpec_mode == "TIMSToF"){
            features <- cbind(features, IM_calibration_vals)
        }

        pth <- base::paste("Temporary_files/Features_aligned_merged",
                           output_file_names_add, ".txt", sep="")
        utils::write.table(features, pth, row.names=FALSE)
        pth <- base::paste("Temporary_files/Feature_alignment_QC_data.RData",
                           sep="")
        save(QC_data, file=pth)
        grDevices::dev.off()
        options(warn=0)
        crap <- gc(FALSE)
    }
}

#' Extend IceR features by expected +1-isotopic features
#' @param path_to_features Path to folder where results of align_features() are
#' stored.
#' @param feature_table_file_name File name which contains align_features()
#' results. By default is set to Features_aligned_merged_IceR_analysis.txt.
#' @param min_observations Specifying how many MaxQuant features had to be
#' previously detected for monoisotopic IceR features to be extended by
#' '+1-isotopic features. By default set to 0 indicating that all for all
#' monoisotopic IceR features a +1-isotopic feature is added.
#' @details Optional step of the IceR workflow appending the list of IceR
#' features by respective expected +1-isotopic features. Isotopic features are
#' expected to show an m/z shift by +1.002 Da/z relative to the monoisotopic
#' IceR feature. Only IceR features with at least min_observations and with mean
#' Andromeda score > 25 % quantile of all mean Andromeda scores are considered.
#' @return Extended list ist stored in the specified feature_table_file_name.
#' @export
add_isotope_features <- function(path_to_features,
                                 feature_table_file_name=def_feat_tab_fname,
                                 min_observations=0){
    setwd(base::paste(path_to_features, "/Temporary_files", sep=""))
    features <- utils::read.table(feature_table_file_name, header=TRUE)
    setwd(path_to_features)

    if(!any(grepl("_i", features$Feature_name))){
        # Isotope + 1 features
        isotope_features <- features
        loc1 <- !grepl(";", isotope_features$Sequence)
        isotope_features <- isotope_features[which(loc1), ]
        loc2 <- as.numeric(as.character(isotope_features$num_matches))
        isotope_features <- isotope_features[which(loc2 >= min_observations), ]

        if(nrow(isotope_features) > 0){
            iso <- as.numeric(as.character(isotope_features$mean_Scores))
            quantile25_score <- grDevices::boxplot.stats(iso)$stats[2]
            loc <- which(iso >= quantile25_score)
            isotope_features <- isotope_features[loc, ]
            isotope_features$m.z <- (((isotope_features$m.z
                                       * isotope_features$Charge) + 1.002054)
                                     / isotope_features$Charge)
            delta_mz <- (isotope_features$m.z_range_max
                         - isotope_features$m.z_range_min)
            dmz <- (delta_mz / 2)
            isotope_features$m.z_range_max <- isotope_features$m.z + dmz
            isotope_features$m.z_range_min <- isotope_features$m.z - dmz
            featn <- base::paste(isotope_features$Feature_name, "_i", sep="")
            isotope_features$Feature_name <- featn
            features <- rbind(features, isotope_features)
            utils::write.table(features, base::paste("Temporary_files/",
                                                     feature_table_file_name,
                                                     sep=""), row.names=FALSE)
            print(paste0(Sys.time(), " Added ", nrow(isotope_features),
                         " isotope features."))
        }

        else{
            print(base::paste(Sys.time(),
                              " None of the peptides was observed in at least",
                              min_observations,
                              "samples. No isotope features were added."))
        }
    }

    else{
        print(paste0(Sys.time(), " Isotope features were already added"))
    }
}

#' Perform quantification of IceR features
#' @param path_to_features Path to folder where results of align_features() are
#' stored
#' @param path_to_mzXML Path to folder containing mzXML files of samples in case
#' of Orbitrap data.
#' @param path_to_MaxQ_output Path to folder containing MaxQuant outputs (txt
#' folder containing at least allpeptides.txt, evidence.txt, peptides.txt and
#' proteinGroups.txt)
#' @param feature_table_file_name File name which contains align_features()
#' results. By default is set to Features_aligned_merged_IceR_analysis.txt.
#' @param output_file_names_add IceR result name tag. By default IceR_analysis
#' @param RT_calibration Boolean value indicating if corrected RT should be
#' used during peak detection, selection and DICE, By default set to T.
#' @param mz_calibration Boolean value indicating if corrected m/z should be
#' used during peak detection, selection and DICE, By default set to T.
#' @param abundance_estimation_correction Boolean value indicating if resulting
#' peptide abundances should be corrected using MaxQuant results as a reference.
#' By default set to T.
#' @param Quant_pVal_cut Numeric value used as diagnostic cutoff border for
# 'visualization of significances of ion accumulation per IceR feature
#' quantification. Furthermore, used as cutoff to filter +1-isotopic IceR
#' features with significant accumulation of ions. By default set to 0.05.
#' @param n_cores Numeric value specifying on how many CPU cores tasks should be
#' parallelized. By default set to 2.
#' @param kde_resolution Numeric value specifying number of grid points per
#' dimension. By default set to 50.
#' @param num_peaks_store Numeric value specifying number of 2D peaks to be
#' stored during peak detection. By default set to 5.
#' @param plot_peak_detection Boolean value indicating if for every feature
#' quantification the determined kernel density estimations and detected peaks
#' should be visualized and stored. By default set to F.
#' @param align_var_score_thr Numeric value specifying
#' significance cutoff to distinguish which features show high RT- or
#' m/z-variability of selected peaks between samples. By default set to 0.05.
#' All features showing significant general variability (variability
#' score < align_var_score_thr) are excluded.
#' @param align_score_thr Numeric value specifying significance cutoff
#' to distinguish which samples show high RT- or m/z-variability of selected
#' peaks for respective IceR feature. By default set to 0.05. All samples
#' showing significant peak variability for respective IceR feature (variability
#' score < align_score_thr) are excluded (quantification set to NA).
#' @param mono_iso_alignment_cutoff Numeric value specifying significance cutoff
#' to distinguish which samples show high RT- or m/z-variability of selected
#' peaks for +1-isotopic from corresponding monoisotopic IceR feature. By
#' default set to 0.05. All samples showing significant peak variability between
#' selected +1-isotopic and monoisotopic IceR features (variability
#' score < mono_iso_alignment_cutoff) are excluded (quantification of
#' '+1-isotopic feature set to NA).
#' @param calc_peptide_LFQ Boolean value specifying if multiply peptide
#' quantification data for same peptide sequence (multiply charge states,
#' isotope-states) should be aggregated using the MaxLFQ algorithm. By default
#' set to F.
#' @param calc_protein_LFQ Boolean value specifying if protein quantification
#' should be additionally performed by peptide quantification aggregation using
#' the MaxLFQ algorithm. By default set to T.
#' @param MassSpec_mode String being either "Orbitrap" or "TIMSToF" specifying
#' by which type of Mass Spectrometer the data was generated. By default it
#' expects Thermo Orbitrap data.
#' @param use_IM_data Boolean value indicating if ion mobility information
#' should be used during feature quantification in case of TIMS-ToF data. By
#' default set to T.
#' @param path_to_extracted_spectra Path to folder containing extracted spectra
#' files of samples in case of TIMS-ToF data.
#' @import foreach
#' @export
#' @details Performs final steps of the IceR workflow:
#' 1) Estimation of background noise per IceR feature quantification.
#' 2) 2D Kernel density estimation-based peak detection, selection and
#' DICE-based quantification of IceR features.
#' 3) Determination of significances of ion accumulations per IceR feature
#' quantification.
#' 4) Quality control of peak selections.
#' 5) IceR peak selection accuracy estimations per sample.
#' 6) Optional: Imputation of missing IceR feature quantifications using
#' estimated backgroudn noise models per sample.
#' 7) Protein quantification by aggregating available peptide quantifications.
#' @return Outputs are stored in the specified output folder and intermediate
#' results are stored in the sub-directory Temporary_files.
#' Quantification results are stored in tab-delimited text files (.tab):
#' Feature information - Features_DDAiceR_Analysis.tab
#' Peak alignment scores - Features_quantification_alignment_score_DDAiceR_
#' Analysis.tab
#' Feature quantifications - Features_quantification_DDAiceR_Analysis.tab
#' Feature quantifications after imputation - Features_quantification_imputed_
#' DDAiceR_Analysis.tab
#' Numbers of observed ions per feature quantification - Features_
#' quantification_ioncount_DDAiceR_Analysis.tab
#' Alignment scores between monoisotopic and corresponding +1isotopic IceR
#' features - Features_quantification_mono_iso_alignment_score_DDAiceR_
#' Analysis.tab
#' Significance of ion accumulations per feature quantification - Features_
#' quantification_pvals_DDAiceR_Analysis.tab
#' Signal to background intensity ratios per feature quantification - Features_
#' quantification_S2B_DDAiceR_Analysis.tab
#' General variability score of peak selections - Features_quantification_
#' variability_score_DDAiceR_Analysis.tab
#' Protein quantification - Proteins_quantification_LFQ_DDAiceR_Analysis.tab
#' Protein quantification after imputation - Proteins_quantification_LFQ_
#' imputed_DDAiceR_Analysis.tab
#' QC results of background noise estimations are visualized in "Decoy feature
#' quantification parameters.pdf".
#' QC results of peak selections are visualized in "Alignment and quantification
#' scores.pdf".
#' The performance of RT- and m/z-alignments over samples is visualized in
#' "Performance of feature alignment.pdf"
#' Estimation of required peptide abundance correction factors is visualized in
#' "Correct feature abundance estimations Signal_Background_intensity.pdf".
#' Intermediate results of the function are stored in RData files.
requantify_features <- function(path_to_features, path_to_mzXML=NA,
                                path_to_MaxQ_output,
                                feature_table_file_name=def_feat_tab_fname,
                                output_file_names_add="IceR_analysis",
                                RT_calibration=TRUE, mz_calibration=TRUE,
                                abundance_estimation_correction=TRUE,
                                Quant_pVal_cut=0.05, n_cores=2,
                                kde_resolution=50, num_peaks_store=5,
                                plot_peak_detection=FALSE,
                                align_var_score_thr=0.05,
                                align_score_thr=0.05,
                                mono_iso_alignment_cutoff=0.05,
                                calc_peptide_LFQ=FALSE, calc_protein_LFQ=TRUE,
                                MassSpec_mode=c("Orbitrap", "TIMSToF"),
                                use_IM_data=TRUE, path_to_extracted_spectra=NA){
    options(warn=-1)
    multiply_intensity_count <- FALSE
    peak_detection <- TRUE

    # Here relevant qc data is stored and finally saved as RData which can be
    # used for re-generating plots
    QC_data <- list()

    mean_background_intensity <- NA
    sd_background_intensity <- NA
    n_background_intensity <- NA

    setwd(base::paste(path_to_features, "/Temporary_files", sep=""))
    features <- utils::read.table(feature_table_file_name,header=TRUE)
    setwd(path_to_features)

    # Check peak widths and increase very small peak widths while reduce extreme
    # long outlier peak widths
    len <- features$RT_length
    stats <- grDevices::boxplot.stats(len)
    # Lower than 25 % quantile
    features$RT_length[which(len < stats$stats[2])] <- stats$stats[2]
    # If RT length > upper whisker --> shorten to upper whisker
    features$RT_length[which(len > stats$stats[5])] <- stats$stats[5]
    features$RT_length[is.na(len)] <- stats$stats[3]

    if(MassSpec_mode == "TIMSToF"){
        len <- features$Inv_K0_length
        stats <- grDevices::boxplot.stats(len)
        # Lower than 25 % quantile
        features$Inv_K0_length[which(len < stats$stats[2])] <- stats$stats[2]
        # If IM length > upper whisker --> shorten to upper whisker
        features$Inv_K0_length[which(len > stats$stats[5])] <- stats$stats[5]
        features$Inv_K0_length[is.na(len)] <- stats$stats[3]
    }

    # Add decoy features
    feats_decoy <- features[which(!grepl("_pmp|_i", features$Feature_name)), ]
    feats_decoy$Feature_name <- base::paste(feats_decoy$Feature_name, "_d",
                                            sep="")
    RT_ranges <- feats_decoy$RT_range_max - feats_decoy$RT_range_min
    mz_ranges <- feats_decoy$m.z_range_max - feats_decoy$m.z_range_min
    median_RT_window <- stats::median(RT_ranges, na.rm=TRUE)
    median_mz_window <- stats::median(mz_ranges, na.rm=TRUE)

    feats_decoy$RT <- feats_decoy$RT + 5 * median_RT_window
    feats_decoy$m.z <- feats_decoy$m.z + 5 * median_mz_window

    feats_decoy$m.z_range_min <- feats_decoy$m.z - (mz_ranges / 2)
    feats_decoy$m.z_range_max <- feats_decoy$m.z + (mz_ranges / 2)

    feats_decoy$RT_range_min <- feats_decoy$RT - (RT_ranges / 2)
    feats_decoy$RT_range_max <- feats_decoy$RT + (RT_ranges / 2)

    # Adjust observed RT and mz accordingly
    Observed_RT <- as.matrix(stringr::str_split(feats_decoy$Observed_RT, ";",
                                                simplify=TRUE))
    class(Observed_RT) <- "numeric"
    Observed_RT <- Observed_RT + 5 * median_RT_window
    feats_decoy$Observed_RT <- apply(Observed_RT, 1, base::paste, collapse=";")

    Observed_mz <- as.matrix(stringr::str_split(feats_decoy$Observed_mz, ";",
                                                simplify=TRUE))
    class(Observed_mz) <- "numeric"
    Observed_mz <- Observed_mz + 5 * median_mz_window
    feats_decoy$Observed_mz <- apply(Observed_mz, 1, base::paste, collapse=";")

    if(MassSpec_mode == "TIMSToF"){
        IM_ranges <- feats_decoy$Inv_K0_range_max-feats_decoy$Inv_K0_range_min
        median_IM_window <- stats::median(IM_ranges, na.rm=TRUE)
        feats_decoy$Inv_K0 <- feats_decoy$Inv_K0 + 5 * median_IM_window

        Observed_IM <- as.matrix(stringr::str_split(feats_decoy$Observed_IM,
                                                    ";", simplify=TRUE))
        class(Observed_IM) <- "numeric"
        Observed_IM <- Observed_IM + 5 * median_IM_window
        feats_decoy$Observed_IM <- apply(Observed_IM, 1, base::paste,
                                         collapse=";")
    }

    features <- rbind(features, feats_decoy)

    if(output_file_names_add != ""){
        output_file_names_add <- base::paste("_", output_file_names_add, sep="")
    }

    # TODO: Move functions to top level to avoid recursive definition
    # Step 1 - Summarize ion intensities per aligned MaxQuant feature
    # Function to perform extraction on multiple threads. Depending on the
    # available ram, the extraction can be run on several threads. However, the
    # task is using much memory so that it is recommended to just run 2 threads
    # in parallel if only 16 gb of ram are available.
    extract_intensities_worker <- function(Sample_IDs, features_select,
                                           path_to_raw, path_to_output_folder,
                                           RT_calibration, mz_calibration,
                                           peak_detection, n_cores,
                                           ion_intensity_cutoff=FALSE,
                                           mean_bkgr_ion_inty_model=NA,
                                           sd_bkgr_ion_inty=NA,
                                           peak_min_ion_count=NA,
                                           kde_resolution=25, num_peaks_store=5,
                                           plots=FALSE,
                                           MassSpec_mode="Orbitrap",
                                           use_IM_data=TRUE){
        # Function to extract intensities in respective extracted .RData for a
        # list of selected features (or all features). indexing_RT_window
        # defines how large each RT indexing window is to speed up subsetting.
        # 0.5 min was observed to be good
        get_intensities <- function(Sample_ID, path, features_select,
                                    indexing_RT_window=0.1, RT_calibration,
                                    mz_calibration, peak_detection,
                                    ion_intensity_cutoff,
                                    mean_bkgr_ion_inty_model,
                                    sd_bkgr_ion_inty,
                                    peak_min_ion_count, kde_resolution,
                                    num_peaks_store, plots, MassSpec_mode,
                                    use_IM_data){
            if(any(!is.na(mean_bkgr_ion_inty_model))){
                cnames <- base::gsub("-", ".",
                                     colnames(mean_bkgr_ion_inty_model))
                colnames(mean_bkgr_ion_inty_model) <- cnames
                rnames <- base::gsub("-", ".",
                                     names(sd_bkgr_ion_inty))
                names(sd_bkgr_ion_inty) <- rnames
            }

            # If no peak selection should be performed, extract all ions within
            # the expected feature windows
            feature_no_peak_selection <- function(all_ion_data,
                                                  selected_features, cur_sample,
                                                  incl_isotope_patterns=FALSE,
                                                  num_peaks_store=0,
                                                  MassSpec_mode="Orbitrap",
                                                  use_IM_data=TRUE){
                peak_ion_data_list <- list()

                # Define RT and mz window
                sel_rt <- selected_features[, base::paste("RT_calibration.",
                                                          cur_sample, sep="")]
                sel_mz <- selected_features[, base::paste("mz_calibration.",
                                                          cur_sample, sep="")]
                RT_expected <- selected_features$RT + sel_rt
                mz_expected <- selected_features$m.z + sel_mz

                RT_window <- c(RT_expected - (selected_features$RT_length / 2),
                               RT_expected + (selected_features$RT_length / 2))

                mz_window <- c(selected_features$m.z_range_min + sel_mz,
                               selected_features$m.z_range_max + sel_mz)

                if(MassSpec_mode == "TIMSToF"){
                    sel <- selected_features[, base::paste("IM_calibration.",
                                                           cur_sample, sep="")]
                    im_expected <- selected_features$Inv_K0 + sel
                    len <- selected_features$Inv_K0_length / 2
                    IM_window <- c(im_expected - len, im_expected + len)
                }

                if(MassSpec_mode == "Orbitrap"){
                    # Select relevant ions
                    sel <- which(all_ion_data$m.z >= mz_window[1] &
                                 all_ion_data$m.z <= mz_window[2] &
                                 all_ion_data$RT >= RT_window[1] &
                                 all_ion_data$RT <= RT_window[2])
                    ion_data <- all_ion_data[sel, ]
                    ion_data <- stats::na.omit(cbind(ion_data,
                                                     data.frame(isotope=0)))

                    # Add isotope ions if wanted
                    if(incl_isotope_patterns == TRUE){
                        iso_mult <- (1:3) * 1.002054
                        cur_charge <- selected_features$Charge
                        cur_mz_min = mz_window[1]
                        cur_mz_max = mz_window[2]

                        for(isos in c(1, 2, 3)){ # Isotope +1,+2,+3
                            c1 <- ((cur_mz_min * cur_charge) + (iso_mult[isos]))
                            c2 <- ((cur_mz_max * cur_charge) + (iso_mult[isos]))
                            mz_range <- c(c1 / cur_charge, c2 / cur_charge)

                            selection <- which(all_ion_data$m.z >= mz_range[1] &
                                               all_ion_data$m.z <= mz_range[2] &
                                               all_ion_data$RT >= RT_window[1] &
                                               all_ion_data$RT <= RT_window[2])

                            if(length(selection) > 0){
                                temp <- all_ion_data[selection, ]
                                temp$isotope <- isos
                                ion_data <- rbind(ion_data, temp)
                            }

                            else{
                                break
                            }
                        }

                        ion_data <- ion_data[order(ion_data$RT), ]
                        # Only keep isotope ions which are present in spectra
                        # where also the previous isotope ions was detected
                        for(iso in c(1, 2, 3)){
                            c1 <- length(which(ion_data$isotope == iso)) > 0
                            c2 <- length(which(ion_data$isotope == iso - 1)) > 0

                            if(c1 & c2){
                                iso_current <- which(ion_data$isotope == iso)

                                loc <- which(ion_data$isotope == iso - 1)
                                iso_minus1 <- ion_data[loc, ]

                                remove <- which(ion_data$RT[iso_current]
                                                %not in% iso_minus1$RT)

                                if(length(remove) > 0){
                                    iso_current <- iso_current[remove]
                                }

                                if(length(iso_current) > 0){
                                    ion_data <- ion_data[-iso_current, ]
                                }
                            }

                            else{
                                remove <- which(ion_data$isotope == iso)

                                if(length(remove) > 0){
                                    ion_data <- ion_data[-remove, ]
                                }
                            }
                        }

                        # Correct m/z of +1,+2 and +3 isotope ions to its
                        # theoretical +0 isotope m/z
                        dif <- ((ion_data$m.z * cur_charge)
                                - (ion_data$isotope * 1.002054))
                        ion_data$m.z <- dif / cur_charge
                    }

                    loc <- as.character(num_peaks_store + 1)
                    df <- base::data.frame(RT=RT_expected,
                                           RT_win_lower=RT_window[1],
                                           RT_win_upper=RT_window[2],
                                           mz=mz_expected,
                                           mz_window_lower=mz_window[1],
                                           mz_window_upper=mz_window[2],
                                           density=0, known_peak=0,
                                           Peak=num_peaks_store + 1,
                                           ion_count=nrow(ion_data))
                    peak_ion_data_list[[loc]] <- list(ion_data=ion_data,
                                                      Peak_info=df)
                }

                if(MassSpec_mode == "TIMSToF"){
                    # Select relevant ions
                    if(use_IM_data == TRUE){
                        loc <- which(all_ion_data$m.z >= mz_window[1] &
                                     all_ion_data$m.z <= mz_window[2] &
                                     all_ion_data$RT >= RT_window[1] &
                                     all_ion_data$RT <= RT_window[2] &
                                     all_ion_data$`1/K0` >= IM_window[1] &
                                     all_ion_data$`1/K0` <= IM_window[2])
                        ion_data <- all_ion_data[loc, ]
                    }

                    else{
                        loc <- which(all_ion_data$m.z >= mz_window[1] &
                                     all_ion_data$m.z <= mz_window[2] &
                                     all_ion_data$RT >= RT_window[1] &
                                     all_ion_data$RT <= RT_window[2])
                        ion_data <- all_ion_data[loc, ]
                    }

                    ion_data <- stats::na.omit(cbind(ion_data,
                                                     data.frame(isotope=0)))

                    loc <- as.character(num_peaks_store + 1)
                    df <- base::data.frame(RT=RT_expected,
                                           RT_win_lower=RT_window[1],
                                           RT_win_upper=RT_window[2],
                                           mz=mz_expected,
                                           mz_window_lower=mz_window[1],
                                           mz_window_upper=mz_window[2],
                                           density=0, known_peak=0,
                                           Peak=num_peaks_store + 1,
                                           ion_count=nrow(ion_data))
                    peak_ion_data_list[[loc]] <- list(ion_data=ion_data,
                                                      Peak_info=df)
                }

                names(peak_ion_data_list) <- "Standard"
                return(peak_ion_data_list)
            }

            # If peak selection should be performed
            feature_2D_peak_selection <- function(all_ion_data,
                                                  selected_features, cur_sample,
                                                  known_RT, delta_mz, delta_rt,
                                                  RT_window_expand_factor=5,
                                                  IM_window_expand_factor=5,
                                                  mz_window_expand_factor=4,
                                                  incl_isotope_patterns=FALSE,
                                                  n_raster_dens_matrix=50,
                                                  local_maxima_k=3,
                                                  max_delta_RT=2,
                                                  max_delta_mz=0.005,
                                                  peak_min_ion_count=5,
                                                  RT_bw=0.5, mz_bw=0.002,
                                                  num_peaks_store=5, plot=FALSE,
                                                  auto_adjust_kde_resol=TRUE,
                                                  MassSpec_mode="Orbitrap",
                                                  delta_im=0.001,
                                                  close_peak_merging=FALSE,
                                                  use_IM_data=TRUE){
                peak_ion_data_list <- list()

                graph <- NA # Here we store graphical output if wanted

                # Define RT and mz window
                n1 <- base::paste("RT_calibration.", cur_sample, sep="")
                n2 <- base::paste("mz_calibration.", cur_sample, sep="")
                RT_correction <- as.numeric(selected_features[n1])
                mz_correction <- as.numeric(selected_features[n2])

                RT_expected <- as.numeric(selected_features$RT + RT_correction)
                mz_expected <- as.numeric(selected_features$m.z + mz_correction)

                # Use standard window as long as the peak width is < standard RT
                # window
                if(delta_rt > selected_features$RT_length / 2){
                    mn <- selected_features$RT_range_min + RT_correction
                    mx <- selected_features$RT_range_max + RT_correction
                    RT_window <- c(mn, mx)
                }

                else{
                    len <- selected_features$RT_length / 2
                    rt <- selected_features$RT + RT_correction
                    RT_window <- c(rt - len, rt + len)
                }

                mz_window <- c(selected_features$m.z_range_min + mz_correction,
                               selected_features$m.z_range_max + mz_correction)

                w1 <- RT_window[1] - ((delta_rt) * RT_window_expand_factor)
                w2 <- RT_window[2] + ((delta_rt) * RT_window_expand_factor)
                RT_window_expanded <- c(w1, w2)

                # Maximum 10 min deviation
                if(RT_window_expanded[2] - RT_window_expanded[1] > 10){
                    RT_window_expanded[1] <- RT_window[1] - 5
                    RT_window_expanded[2] <- RT_window[2] + 5
                }

                xpnd <- (delta_mz * mz_window_expand_factor)
                mz_window_expanded <- c(mz_window[1] - xpnd,
                                        mz_window[2] + xpnd)

                if(MassSpec_mode == "TIMSToF"){
                    loc <- base::paste("IM_calibration.", cur_sample, sep="")
                    IM_correction <- as.numeric(selected_features[loc])
                    IM_expected <- as.numeric(selected_features$Inv_K0
                                              + IM_correction)

                    # Use standard window as long as the peak width
                    # is < standard IM window
                    if(delta_im > selected_features$Inv_K0_length / 2){
                        mn <- selected_features$Inv_K0_range_min + IM_correction
                        mx <- selected_features$Inv_K0_range_max + IM_correction
                        IM_window <- c(mn, mx)
                    }

                    else{
                        k0_cor <- selected_features$Inv_K0 + IM_correction
                        len <- (selected_features$Inv_K0_length / 2)
                        IM_window <- c(k0_cor - len, k0_cor + len)
                    }

                    dt <- (delta_im * IM_window_expand_factor)
                    IM_window_expanded <- c(IM_window[1] - dt,
                                            IM_window[2] + dt)
                }

                # Filter for ions in the expanded window
                if(MassSpec_mode == "Orbitrap"){
                    loc <- which(all_ion_data$m.z >= mz_window_expanded[1] &
                                 all_ion_data$m.z <= mz_window_expanded[2] &
                                 all_ion_data$RT >= RT_window_expanded[1] &
                                 all_ion_data$RT <= RT_window_expanded[2])
                    ion_data <- all_ion_data[loc, ]
                }

                else{
                    if(use_IM_data == TRUE){
                        loc <- which(all_ion_data$m.z >= mz_window_expanded[1] &
                                     all_ion_data$m.z <= mz_window_expanded[2] &
                                     all_ion_data$RT >= RT_window_expanded[1] &
                                     all_ion_data$RT <= RT_window_expanded[2] &
                                     all_ion_data$`1/K0` >= IM_window[1] &
                                     all_ion_data$`1/K0` <= IM_window[2])
                        ion_data <- all_ion_data[loc, ]
                    }

                    else{
                        loc <- which(all_ion_data$m.z >= mz_window_expanded[1] &
                                     all_ion_data$m.z <= mz_window_expanded[2] &
                                     all_ion_data$RT >= RT_window_expanded[1] &
                                     all_ion_data$RT <= RT_window_expanded[2])
                        ion_data <- all_ion_data[loc, ]
                    }
                }

                ion_data <- stats::na.omit(cbind(ion_data,
                                                 data.frame(isotope=0)))

                # Add isotope ions if wanted
                if(incl_isotope_patterns == TRUE & nrow(ion_data) > 0){
                    iso_mult <- (1:3) * 1.002054

                    cur_charge <- selected_features$Charge
                    cur_mz_min <- mz_window_expanded[1]
                    cur_mz_max <- mz_window_expanded[2]

                    # Isotope +1,+2,+3
                    for(isos in c(1, 2, 3)){
                        mz_range <- c(((cur_mz_min * cur_charge)
                                       + (iso_mult[isos])) / cur_charge,
                                      ((cur_mz_max * cur_charge)
                                       + (iso_mult[isos])) / cur_charge)

                        sel <- which(all_ion_data$m.z >= mz_range[1] &
                                     all_ion_data$m.z <= mz_range[2] &
                                     all_ion_data$RT >= RT_window_expanded[1] &
                                     all_ion_data$RT <= RT_window_expanded[2])

                        if(length(sel) > 0){
                            temp <- all_ion_data[sel, ]
                            temp$isotope <- isos
                            ion_data <- rbind(ion_data, temp)
                        }

                        else{
                            break
                        }
                    }

                    ion_data <- ion_data[order(ion_data$RT), ]
                    # Only keep isotope ions which are present in spectra where
                    # also the previous isotope ions was detected
                    for(iso in c(1, 2, 3)){
                        c1 <- length(which(ion_data$isotope == iso)) > 0
                        c2 <- length(which(ion_data$isotope == iso - 1)) > 0

                        if(c1 & c2){
                            iso_current <- which(ion_data$isotope == iso)
                            loc <- which(ion_data$isotope == iso - 1)
                            iso_minus1 <- ion_data[loc, ]

                            remove <- which(ion_data$RT[iso_current]
                                            %not in% iso_minus1$RT)

                            if(length(remove) > 0){
                                iso_current <- iso_current[remove]
                            }

                            if(length(iso_current) > 0){
                                ion_data <- ion_data[-iso_current, ]
                            }
                        }

                        else{
                            remove <- which(ion_data$isotope == iso)

                            if(length(remove) > 0){
                                ion_data <- ion_data[-remove, ]
                            }
                        }
                    }

                    # Correct m/z of +1,+2 and +3 isotope ions to its
                    # theoretical +0 isotope m/z
                    num <- ((ion_data$m.z * cur_charge)
                            - (ion_data$isotope * 1.002054))
                    ion_data$m.z <- num / cur_charge
                }

                if(nrow(ion_data) >= peak_min_ion_count){
                    # Determine 2D density matrix
                    # Consider different resolutions required for peaks with
                    # expected small RT lengths
                    if(auto_adjust_kde_resol == TRUE){
                        num <- (RT_window_expanded[2] - RT_window_expanded[1])
                        den <- ((selected_features$RT_length / 1.9) * 2)
                        required_res <- ceiling((num / den) + 1)

                        if(required_res > n_raster_dens_matrix){
                            n_raster_dens_matrix <- required_res
                        }
                    }

                    f2 <- MASS::kde2d(ion_data$RT, ion_data$m.z,
                                      n=n_raster_dens_matrix, h=c(RT_bw, mz_bw),
                                      lims=c(RT_window_expanded[1],
                                             RT_window_expanded[2],
                                             mz_window_expanded[1],
                                             mz_window_expanded[2]))

                    # Determine which density (=z) will be selected to be at
                    # least exceeded --> upper whisker
                    outlier_densities <- grDevices::boxplot.stats(f2$z)$stats[5]

                    # Convert it to a raster object
                    r <- raster::raster(f2$z)
                    raster::extent(r) <- raster::extent(c(0, length(f2$x), 0,
                                                          length(f2$y)) + 0.5)

                    # Find the maximum value within the k-cell neighborhood of
                    # each cell
                    f <- function(X){
                        max(X, na.rm=TRUE)
                    }

                    # Weight matrix for cells in moving window
                    ww <- matrix(1, nrow=local_maxima_k, ncol=local_maxima_k)
                    localmax <- raster::focal(r, fun=f, w=ww, pad=TRUE,
                                              padValue=NA)

                    # Get x-y coordinates of those cells that are local maxima
                    r2 <- r == localmax
                    maxXY <- raster::xyFromCell(r2, raster::Which(r2==1,
                                                                  cells=TRUE))

                    # Remove maxima which are below the minimal outlier density
                    rt <- f2$x[n_raster_dens_matrix - maxXY[, 2] + 1]
                    mz <- f2$y[maxXY[, 1]]
                    dty <- f2$z[((maxXY[, 1] - 1) * n_raster_dens_matrix)
                                + (n_raster_dens_matrix - maxXY[, 2] + 1)]
                    maxima <- base::data.frame(RT=rt, mz=mz, density=dty)
                    loc <- which(maxima$density > outlier_densities)
                    maxima <- maxima[loc, ]
                    # Next filter for maxima with at least peak_min_ion_count
                    # ions
                    div <- (maxima$density / sum(as.matrix(maxima$density)))
                    maxima$estimated_count <- div * nrow(ion_data)

                    RT_cut <- selected_features$RT_length / 2
                    c1 <- maxima$estimated_count >= peak_min_ion_count
                    c2 <- abs(maxima$RT - RT_expected) <= RT_cut
                    c3 <- abs(maxima$mz - mz_expected) <= delta_mz
                    selection <- which(c1 | c2 & c3)
                    # If no maxima are left after filtering still take topN
                    # closest peak further
                    if(length(selection) == 0 & nrow(maxima) > 0){
                        # No maxima would be left see keep up to N peaks
                        drt <- maxima$RT - RT_expected
                        dmz <- maxima$mz - mz_expected
                        dts <- sqrt(((drt))^2 + ((dmz) * 500)^2)
                        distances <- base::data.frame(dist_rt=drt, dist_mz=dmz,
                                                      dist_total_scaled=dts)
                        # At maximum the top N closest maxima
                        ordering <- order(abs(distances$dist_total_scaled))
                        maxima <- maxima[ordering, ]

                        if(nrow(maxima) > num_peaks_store){
                            maxima <- maxima[1:num_peaks_store, ]
                        }
                    }

                    else if(nrow(maxima) > 0){
                        maxima <- maxima[selection, ]
                    }

                    # Next merge peaks which are very close together
                    if(nrow(maxima)>1 & close_peak_merging == TRUE){
                        mat <- matrix(ncol=2, nrow=nrow(maxima), 0)
                        close_info <- base::as.data.frame(mat)

                        for(cl in 1:nrow(maxima)){
                            n1 <- (maxima$RT[cl] - maxima$RT)^2
                            n2 <- ((maxima$mz[cl] - maxima$mz) * 500)^2
                            temp <- sqrt(n1 + n2)
                            sel <- which(temp == min(temp[-cl]))[1]
                            data.table::set(close_info, as.integer(cl),
                                            as.integer(1:2),
                                            value=list(temp[sel], sel))
                        }

                        # Merge peaks with d_mz < 0.001 and d_RT < peak_width/2
                        n1 <- (delta_mz * 500)^2
                        n2 <- (selected_features$RT_length / 4)^2
                        dist_cut <- sqrt(n1 + n2)

                        if(any(close_info$V1 < dist_cut)){
                            sel <- which(close_info$V1 < dist_cut)

                            mat <- matrix(nrow=length(sel), ncol=4, 0)
                            maxima_add <- base::as.data.frame(mat)
                            colnames(maxima_add) <- colnames(maxima)
                            # Find mean RT and mz of maxima which should be
                            # merged wheigted by density
                            for(cl in 1:length(sel)){
                                ssel <- close_info[sel[cl], 2]
                                rt_s1 <- maxima$RT[sel[cl]]
                                rt_s2 <- maxima$RT[ssel]
                                dt_s1 <- maxima$density[sel[cl]]
                                dt_s2 <- maxima$density[ssel]
                                mz_s1 <- maxima$mz[sel[cl]]
                                mz_s2 <- maxima$mz[ssel]
                                ct_s1 <- maxima$estimated_count[sel[cl]]
                                ct_s2 <- maxima$estimated_count[ssel]
                                RT_add <- stats::weighted.mean(c(rt_s1, rt_s2),
                                                               c(dt_s1, dt_s2))
                                mz_add <- stats::weighted.mean(c(mz_s1, mz_s2),
                                                               c(dt_s1, dt_s2))
                                dty_add <- stats::weighted.mean(c(dt_s1, dt_s2),
                                                                c(dt_s1, dt_s2))
                                cnt_add <- stats::weighted.mean(c(ct_s1, ct_s2),
                                                                c(dt_s1, dt_s2))

                                data.table::set(maxima_add, as.integer(cl),
                                                as.integer(1:4),
                                                value=list(RT_add,　mz_add,
                                                        　　dty_add, cnt_add))
                            }
                            maxima_add <- unique(maxima_add)
                            # Remove maxima which are merged
                            remove <- unique(append(sel, close_info[sel, 2]))
                            maxima <- maxima[-remove, ]
                            # Add new merged peaks
                            maxima <- rbind(maxima, maxima_add)
                        }
                    }

                    #　Now check which density maxima is closest to the expected
                    # feature
                    if(nrow(maxima) > 0){
                        drt <- maxima$RT - RT_expected
                        dmz <- maxima$mz - mz_expected
                        dts <- sqrt(((maxima$RT - RT_expected))^2
                                    + ((maxima$mz - mz_expected) * 500)^2)
                        distances <- base::data.frame(dist_rt=drt, dist_mz=dmz,
                                                      dist_total_scaled=dts)
                        # At maximum the top N closest maxima
                        ordering <- order(abs(distances$dist_total_scaled))
                        maxima_select <- maxima[ordering, ]

                        if(nrow(maxima_select) > num_peaks_store){
                            maxima_select <- maxima_select[1:num_peaks_store,　]
                        }

                        if(nrow(maxima_select)　>　0){
                            # Label which peak was closest to expected feature
                            # RT and mz if known
                            maxima_select$known_peak <- 0

                            if(!is.na(known_RT)){
                                # RT is known
                                # Check if any peak was detected close to known
                                # RT (max deviation by peak width/2)
                                dif1 <- abs(maxima_select$RT - RT_expected)
                                dif2 <- abs(maxima_select$mz - mz_expected)
                                c1 <- dif1 <= selected_features$RT_length / 1.9
                                c2 <- dif2 <= delta_mz / 1.9
                                peak_close_rt <- which(c1 & c2)
                                # And if the closest peak actually also shows a
                                # good number of ions. Otherwise no peak will be
                                # selected here but selection will be done later

                                c1 <- length(peak_close_rt) > 0
                                c2 <- (maxima_select$estimated_count[1] >=
                                       2 * peak_min_ion_count)

                                if(c1 & c2){
                                    pk <- peak_close_rt[1]
                                    maxima_select$known_peak[pk] <- 1
                                    # Reorder maxima selected such that known
                                    # peak maximum is at index 1
                                    loc <- c(pk, c(1:nrow(maxima_select))[-pk])
                                    maxima_select <- maxima_select[loc, ]
                                }
                                # Peak RT is known but none of the detected
                                # peaks shows RT close to the known RT
                                else{
                                    # In this case none of the maxima are set to
                                    # be the known peak and later we will try to
                                    # decide which one is the correct one
                                }
                            }
                            # If any maximum was detected, return table of ions
                            # as a list containing up to 3 peak area ion lists
                            # from closest to further away peaks.
                            RT_window_width <- selected_features$RT_length
                            mz_window_width <- delta_mz

                            peak_counter <- 0

                            for(i in 1:nrow(maxima_select)){
                                mzi <- maxima_select$mz[i]
                                rti <- maxima_select$RT[i]
                                cur_mz_window <- c(mzi - (mz_window_width),
                                                   mzi + (mz_window_width))
                                cur_RT_win <- c(rti - (RT_window_width / 2),
                                                   rti + (RT_window_width / 2))

                                mz <- ion_data$m.z
                                rt <- ion_data$RT
                                selection <- which(mz >= cur_mz_window[1] &
                                                   mz <= cur_mz_window[2] &
                                                   rt >= cur_RT_win[1] &
                                                   rt <= cur_RT_win[2])

                                c1 <- length(selection) < peak_min_ion_count
                                c2 <- maxima_select$known_peak[i] == 1

                                if(c1 & c2){
                                    # If this peak should be the correct peak
                                    # but number of ions are very low then set
                                    # to unknown correct peak
                                    maxima_select$known_peak[i] <- 0
                                }

                                peak_counter <- peak_counter + 1
                                df <- base::data.frame(
                                    RT=maxima_select$RT[i],
                                    RT_win_lower=cur_RT_win[1],
                                    RT_win_upper=cur_RT_win[2],
                                    mz=maxima_select$mz[i],
                                    mz_window_lower=cur_mz_window[1],
                                    mz_window_upper=cur_mz_window[2],
                                    density=maxima_select$density[i],
                                    known_peak=maxima_select$known_peak[i],
                                    Peak=peak_counter,
                                    ion_count=length(selection)
                                )
                                lst <- list(ion_data=ion_data[selection, ],
                                            Peak_info=df)
                                peak_ion_data_list[[peak_counter]] <- lst
                            }

                            if(peak_counter > 0){
                                nms <- as.character(1:peak_counter)
                                names(peak_ion_data_list) <- nms
                            }

                            w <- (RT_window_width / 2)
                            cur_RT_win <- c(RT_expected - w,
                                               RT_expected + w)
                            # Also add quantification if standard windows are
                            # used
                            selection <- which(ion_data$m.z >= mz_window[1] &
                                               ion_data$m.z <= mz_window[2] &
                                               ion_data$RT >= cur_RT_win[1] &
                                               ion_data$RT <= cur_RT_win[2])
                            loc <- as.character(num_peaks_store + 1)
                            df <- base::data.frame(RT=RT_expected,
                                                   RT_win_lower=cur_RT_win[1],
                                                   RT_win_upper=cur_RT_win[2],
                                                   mz=mz_expected,
                                                   mz_window_lower=mz_window[1],
                                                   mz_window_upper=mz_window[2],
                                                   density=0, known_peak=0,
                                                   Peak=num_peaks_store + 1,
                                                   ion_count=length(selection))
                            lst <- list(ion_data=ion_data[selection, ],
                                        Peak_info=df)
                            peak_ion_data_list[[loc]] <- lst
                            # Plot should be directly generated
                            if(plot == TRUE){
                                feats <- selected_features$Feature_name
                                title <- base::paste(cur_sample, "-", feats)
                                graphics::image(f2, main=title, xlab="RT",
                                                ylab="m/z",
                                                xlim=RT_window_expanded,
                                                ylim=mz_window_expanded)
                                graphics::text(maxima_select$RT,
                                               maxima_select$mz,
                                               1:nrow(maxima_select),
                                               col="black")
                                # Expected window
                                graphics::rect(cur_RT_win[1], mz_window[1],
                                               cur_RT_win[2], mz_window[2],
                                               lty=2)
                                graphics::points(RT_expected, mz_expected,
                                                 pch=4)
                                # Adjusted window
                                dif1 <- maxima_select$RT - RT_window_width / 2
                                dif2 <- maxima_select$mz - mz_window_width
                                sum1 <- maxima_select$RT + RT_window_width / 2
                                sum2 <- maxima_select$mz + mz_window_width
                                border <- ifelse(maxima_select$known_peak == 1,
                                                 "green", "darkgrey")
                                graphics::rect(dif1, dif2, sum1, sum2, lty=2,
                                               border=border)

                                # Add label indicating intensity within the
                                # window
                                for(i in 1:length(peak_ion_data_list)){
                                    pki <- peak_ion_data_list[[i]]
                                    sum_inty <- sum(sum(pki$ion_data$Intensity))

                                    if(sum_inty > 0){
                                        sum_inty <- base::log2(sum_inty)
                                    }

                                    pkinf <- pki$Peak_info
                                    dt <- (1.1 * delta_mz)
                                    col <- ifelse(pkinf$known_peak == 1,
                                                  "green", "darkgrey")

                                    if(pkinf$Peak != num_peaks_store + 1){
                                        graphics::text(pkinf$RT, pkinf$mz - dt,
                                                       round(sum_inty, 2),
                                                       col=col)
                                    }

                                    else{
                                        graphics::text(pkinf$RT, pkinf$mz + dt,
                                                       round(sum_inty, 2),
                                                       col="black")
                                    }
                                }
                            }

                            # Store all relevant data for performing plotting
                            # later
                            else{
                                mat <- matrix(ncol=1,
                                              nrow=length(peak_ion_data_list))
                                peak_intensities <- base::as.data.frame(mat)
                                pk <- as.numeric(peak_intensities[, 1])
                                peak_intensities[, 1] <- pk

                                for(i in 1:length(peak_ion_data_list)){
                                    pki <- peak_ion_data_list[[i]]$ion_data
                                    sum_inty <- sum(sum(pki$Intensity))

                                    if(sum_inty > 0){
                                        sum_inty <- base::log2(sum_inty)
                                    }

                                    data.table::set(peak_intensities,
                                                    as.integer(i),
                                                    as.integer(1), sum_inty)
                                }

                                graph <- list(kdemap=f2,
                                              maxima_select=maxima_select,
                                              peak_intensities=peak_intensities)
                            }
                        }

                        # No maxima with a minmal density found
                        else{
                            #In this case simply extract ions which are within
                            # the expected window
                            RT_window_width <- selected_features$RT_length
                            rw <- (RT_window_width / 2)
                            cur_RT_win <- c(RT_expected - rw, RT_expected + rw)

                            selection <- which(ion_data$m.z >= mz_window[1] &
                                               ion_data$m.z <= mz_window[2] &
                                               ion_data$RT >= cur_RT_win[1] &
                                               ion_data$RT <= cur_RT_win[2])
                            loc <- as.character(num_peaks_store + 1)
                            df <- base::data.frame(RT=RT_expected,
                                                   RT_win_lower=cur_RT_win[1],
                                                   RT_win_upper=cur_RT_win[2],
                                                   mz=mz_expected,
                                                   mz_window_lower=mz_window[1],
                                                   mz_window_upper=mz_window[2],
                                                   density=0, known_peak=0,
                                                   Peak=4,
                                                   ion_count=length(selection))
                            lst <- list(ion_data=ion_data[selection, ],
                                        Peak_info=df)
                            peak_ion_data_list[[loc]] <- lst
                        }
                    }

                    # No maxima found
                    else{
                        # In this case simply extract ions which are within the
                        # expected window
                        RT_window_width <- selected_features$RT_length
                        rw <- (RT_window_width / 2)
                        cur_RT_win <- c(RT_expected - rw, RT_expected + rw)

                        selection <- which(ion_data$m.z >= mz_window[1] &
                                           ion_data$m.z <= mz_window[2] &
                                           ion_data$RT >= cur_RT_win[1] &
                                           ion_data$RT <= cur_RT_win[2])
                        loc <- as.character(num_peaks_store + 1)
                        df <- base::data.frame(RT=RT_expected,
                                               RT_win_lower=cur_RT_win[1],
                                               RT_win_upper=cur_RT_win[2],
                                               mz=mz_expected,
                                               mz_window_lower=mz_window[1],
                                               mz_window_upper=mz_window[2],
                                               density=0, known_peak=0,
                                               Peak=4,
                                               ion_count=length(selection))
                        lst <- list(ion_data=ion_data[selection, ],
                                    Peak_info=df)
                        peak_ion_data_list[[loc]] <- lst
                    }
                }

                # Not enough ions even in expanded window
                else{
                    # In this case simply extract ions which are within the
                    # expected window
                    RT_window_width <- selected_features$RT_length
                    cur_RT_win <- c(RT_expected - (RT_window_width / 2),
                                    RT_expected + (RT_window_width / 2))

                    selection <- which(ion_data$m.z >= mz_window[1] &
                                       ion_data$m.z <= mz_window[2] &
                                       ion_data$RT >= cur_RT_win[1] &
                                       ion_data$RT <= cur_RT_win[2])
                    loc <- as.character(num_peaks_store + 1)
                    df <- base::data.frame(RT=RT_expected,
                                           RT_win_lower=cur_RT_win[1],
                                           RT_win_upper=cur_RT_win[2],
                                           mz=mz_expected,
                                           mz_window_lower=mz_window[1],
                                           mz_window_upper=mz_window[2],
                                           density=0, known_peak=0, Peak=4,
                                           ion_count=length(selection))
                    lst <- list(ion_data=ion_data[selection, ],
                                Peak_info=df)
                    peak_ion_data_list[[loc]] <- lst
                }

                npks <- as.character(num_peaks_store + 1)
                nms <- names(peak_ion_data_list)
                loc1 <- which(nms != npks)
                loc2 <- which(nms == npks)
                names(peak_ion_data_list)[loc1] <- base::paste("Peak_",
                                                                nms[loc1],
                                                                sep="")
                names(peak_ion_data_list)[loc2] <- "Standard"
                peak_ion_data_list$graph <- graph

                return(peak_ion_data_list)
            }

            `%not in%` <- function(x, table){
                is.na(match(x, table, nomatch=NA_integer_))
            }

            t.test2 <- function(...){
                obj <- try(stats::t.test(...), silent=TRUE)

                if(methods::is(obj, "try-error")){
                    return(NA)
                }

                else{
                    return(obj$p.value)
                }
            }

            setwd(path)

            max <- 1
            print("Load spectra data")
            title <- base::paste("Load spectra data:", Sample_ID)
            label <- base::paste(round(0 / max * 100, 0), "% done")
            pb <- tcltk::tkProgressBar(title=title, label=label, min=0, max=max,
                                       width=300)
            table_store <- NULL

            # SILAC mode --- all channels in one raw file
            grp <- "_Channel_light|_Channel_medium|_Channel_heavy"

            if(grepl(grp, Sample_ID)){
                # Ions with RT and intensity in variable dat
                load(base::paste(base::gsub(grp, "", Sample_ID),
                                 "_all_ions.RData", sep=""))
            }

            # Ions with RT and intensity in variable dat
            else{
                load(base::paste(Sample_ID, "_all_ions.RData", sep=""))
            }

            label <- base::paste(round(1 / max * 100, 0), " % done (", 1, "/",
                                 max, ")", sep="")
            tcltk::setTkProgressBar(pb, 1, label=label)
            close(pb)

            if(MassSpec_mode == "TIMSToF"){
                dat <- table_store
                rm(table_store)
                gc()

                dat$RT <- dat$RT / 60
                dat <- dat[, c(3, 1, 4, 2)]
                colnames(dat) <- c("m.z", "RT", "Intensity", "1/K0")
            }

            # Indexing dat by RT windows to improve speed for subsetting
            cl <- ceiling(max(dat$RT))
            num_windows <- ceiling(cl * (1 / indexing_RT_window))

            # If in TIMSToF mode, further index by IM to further improve speed
            if(MassSpec_mode == "TIMSToF"){
                indexing_IM_window <- 0.025
                cl <- ceiling(max(dat$`1/K0`) - min(dat$`1/K0`))
                num_IM_windows <- ceiling(cl * (1 / indexing_IM_window))
                num_windows_RT_only <- num_windows
                min_im <- min(dat$`1/K0`)
                num_windows <- num_windows * num_IM_windows
            }

            if(MassSpec_mode == "Orbitrap"){
                Indices <- base::as.data.frame(matrix(ncol=4, nrow=num_windows))
                colnames(Indices) <- c("RT_start", "RT_end", "Row_start",
                                       "Row_end")
                Indices$RT_start <- as.numeric(Indices$RT_start)
                Indices$RT_end <- as.numeric(Indices$RT_end)
                Indices$Row_start <- as.numeric(Indices$Row_start)
                Indices$Row_end <- as.numeric(Indices$Row_end)

                print("Indexing intensities")
                max <- nrow(Indices)
                title <- base::paste("Indexing intensities:", Sample_ID)
                label <- base::paste(round(0 / max * 100, 0), "% done")
                pb <- tcltk::tkProgressBar(title=title, label=label, min=0,
                                           max=max, width=300)
                start_time <- Sys.time()
                updatecounter <- 0
                time_require <- 0
                start <- 1
                end <- max
                cur_index <- 1

                # Increment per indexing
                incr <- 10000

                for(i in 1:nrow(Indices)){
                    start <- (i - 1) * indexing_RT_window
                    end <- i * indexing_RT_window

                    indx <- cur_index
                    indx_prev <- cur_index

                    while(TRUE){
                        indx <- indx + incr

                        if(dat$RT[indx] >= end | indx > nrow(dat)){
                            break
                        }

                        else{
                            indx_prev <- indx
                        }
                    }

                    if(indx > nrow(dat)){
                        indx <- nrow(dat)
                    }


                    if(indx > indx_prev){
                        temp <- dat[indx_prev:indx,]

                        if(length(which(temp$RT < end)) > 0){
                            mx <- max(which(temp$RT < end))
                            inds <- cur_index:(indx_prev - 1 + mx)
                            incr <- max(inds) + 1 - cur_index
                            cur_index <- max(inds) + 1

                            if(length(inds) > 0){
                                data.table::set(Indices, i=as.integer(i),
                                                j=as.integer(1:4),
                                                value=as.list(c(start, end,
                                                                min(inds),
                                                                max(inds))))
                            }
                        }
                    }

                    updatecounter <- updatecounter + 1

                    if(updatecounter >= 1){
                        time_elapsed <- difftime(Sys.time(), start_time,
                                                 units="secs")
                        imx <- (i / max)
                        time_require <- (time_elapsed / imx) * (1 - imx)
                        td <- lubridate::seconds_to_period(time_require)
                        time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                                lubridate::minute(td),
                                                round(lubridate::second(td),
                                                digits=0))
                        updatecounter <- 0
                        label <- base::paste(round(imx * 100, 0), " % done (",
                                             i, "/", max, ", Time require: ",
                                             time_require, ")", sep="")
                        tcltk::setTkProgressBar(pb, i, label=label)
                    }

                }

                close(pb)

                if(any(rowSums(is.na(Indices)) > 0)){
                    Indices <- Indices[-which(rowSums(is.na(Indices)) > 0), ]
                }

                mn <- min(features_select$RT_range_min)

                if(mn < Indices$RT_start[1]){
                    Indices$RT_start[1] <- mn
                }

                mx <- max(features_select$RT_range_max)

                if(mx > Indices$RT_end[nrow(Indices)]){
                    Indices$RT_end[nrow(Indices)] <- mx
                }

                # Fragment data into indexed RT windows to further improve
                # subsetting speed
                data_frags <- list()

                for(i in 1:nrow(Indices)){
                    ids <- Indices$Row_start[i]:Indices$Row_end[i]
                    data_frags[[i]] <- dat[ids, ]
                }

                rm(dat)
                crap <- gc(FALSE)
            }

            else{
                Indices <- base::as.data.frame(matrix(ncol=4, nrow=num_windows))
                colnames(Indices) <- c("RT_start", "RT_end", "IM_start",
                                       "IM_end")
                r1 <- (0:(num_windows_RT_only - 1))
                r2 <- (1:num_windows_RT_only)
                Indices$RT_start <- sort(rep((r1 * indexing_RT_window),
                                             num_IM_windows))
                Indices$RT_end <- sort(rep((r2 * indexing_RT_window),
                                           num_IM_windows))
                Indices$IM_start <- min_im + (r1 * indexing_IM_window)
                Indices$IM_end <- min_im + (r2 * indexing_IM_window)

                Indices_list <- list()
                Indices_list_count <- 0

                print("Indexing intensities")
                max <- num_windows_RT_only
                title <- base::paste("Indexing intensities:", Sample_ID)
                label <- base::paste(round(0 / max * 100, 0), "% done")
                pb <- tcltk::tkProgressBar(title=title, label=label, min=0,
                                           max=max, width=300)
                start_time <- Sys.time()
                updatecounter <- 0
                time_require <- 0
                start <- 1
                end <- max

                cur_index <- 1

                # Increment per indexing
                incr <- 10000

                for(i in 1:num_windows_RT_only){
                    start <- (i - 1) * indexing_RT_window
                    end <- i * indexing_RT_window

                    indx <- cur_index
                    indx_prev <- cur_index

                    while(TRUE){
                        indx <- indx + incr

                        if(dat$RT[indx] >= end | indx > nrow(dat)){
                            break
                        }

                        else{
                            indx_prev <- indx
                        }
                    }

                    if(indx > nrow(dat)){
                        indx <- nrow(dat)
                    }

                    if(indx > indx_prev){
                        temp <- dat[indx_prev:indx, ] # XXX: FROM HERE

                        if(length(which(temp$RT < end)) > 0){
                            mx <- max(which(temp$RT < end))
                            inds <- cur_index:(indx_prev - 1 + mx)
                            incr <- max(inds) + 1 - cur_index
                            temp <- dat[inds, ]

                            # Now further subset based on IM
                            for(j in 1:num_IM_windows){
                                start_im <- min_im + ((j - 1)
                                                      * indexing_IM_window)
                                end_im <- min_im + (j * indexing_IM_window)
                                loc <- which(temp$`1/K0` >= start_im
                                             & temp$`1/K0` <= end_im)
                                itemp <- inds[loc]
                                Indices_list_count <- Indices_list_count + 1

                                if(length(itemp) > 0){
                                    Indices_list[[Indices_list_count]] <- itemp
                                }

                                else{
                                    Indices_list[[Indices_list_count]] <- NA
                                }
                            }

                            cur_index <- max(inds) + 1
                        }
                    }

                    if(Indices_list_count != (i * num_IM_windows)){
                        missing <- (i * num_IM_windows) - Indices_list_count

                        for(j in 1:missing){
                            Indices_list_count <- Indices_list_count + 1
                            Indices_list[[Indices_list_count]] <- NA
                        }
                    }

                    updatecounter <- updatecounter + 1

                    if(updatecounter >= 1){
                        time_elapsed <- difftime(Sys.time(), start_time,
                                                 units="secs")
                        time_require <- ((time_elapsed / (i / max))
                                         * (1 - (i / max)))
                        td <- lubridate::seconds_to_period(time_require)
                        time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                                lubridate::minute(td),
                                                round(lubridate::second(td),
                                                      digits=0))
                        updatecounter <- 0
                        label <- base::paste(round(i / max * 100, 0),
                                             " % done (", i, "/", max,
                                             ", Time require: ", time_require,
                                             ")", sep = "")
                        tcltk::setTkProgressBar(pb, i, label=label)
                    }
                }

                close(pb)

                if(any(is.na(Indices_list))){
                    remove <- which(is.na(Indices_list))
                    Indices <- Indices[-remove, ]

                    na.omit.list <- function(y){
                        return(y[!sapply(y, function(x) all(is.na(x)))])
                    }

                    Indices_list <- na.omit.list(Indices_list)
                }

                # Fragment data into indexed RT windows to further improve
                # subsetting speed
                data_frags <- list()

                print("Fragment data")
                max <- nrow(Indices)
                title <- base::paste("Fragment data:", Sample_ID)
                label <- base::paste(round(0 / max * 100, 0), "% done")
                pb <- tcltk::tkProgressBar(title=title, label=label, min=0,
                                           max=max, width=300)
                start_time <- Sys.time()
                updatecounter <- 0
                time_require <- 0

                for(i in 1:nrow(Indices)){
                    data_frags[[i]] <- dat[Indices_list[[i]], ]
                    updatecounter <- updatecounter + 1

                    if(updatecounter >= 100){
                        time_elapsed <- difftime(Sys.time(), start_time,
                                                 units="secs")
                        time_require <- ((time_elapsed / (i / max))
                                         * (1 - (i / max)))
                        td <- lubridate::seconds_to_period(time_require)
                        time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                                lubridate::minute(td),
                                                round(lubridate::second(td),
                                                      digits=0))

                        updatecounter <- 0
                        label <- base::paste(round(i / max * 100, 0),
                                             " % done (", i, "/", max,
                                             ", Time require: ", time_require,
                                             ")", sep="")
                        tcltk::setTkProgressBar(pb, i, label=label)
                    }
                }

                close(pb)
            }

            rm(dat)
            crap <- gc(FALSE)

            # Now extract intensities,
            # row1 = Intensity summed,
            # row2 = num ions,
            # row3 = mean intensity,
            # row4 = sd intensity,
            # row5 = detected optimum mz window minimum
            # row6 = detected optimum mz window maximum
            # row7 = detected optimum RT window minimum
            # row8 = detected optimum RT window maximum
            # row9 = number of detected peaks

            if(peak_detection == FALSE){
                graph_peaks <- list()
                peaks_quant <- list()

                mat <- matrix(nrow=10, ncol=nrow(features_select), 0)
                Intensities <- base::as.data.frame(mat)
                colnames(Intensities) <- features_select$Feature_name

                # If cut off pvalue is defined then store background
                # quantifications separately
                if(ion_intensity_cutoff == TRUE){
                    mat <- matrix(nrow=10, ncol=nrow(features_select), 0)
                    Intensities_signal_background <- base::as.data.frame(mat)
                    ftn <- features_select$Feature_name
                    colnames(Intensities_signal_background) <- ftn
                }

                # Store data for calculating scoring per sample and feature
                mat <- matrix(nrow=1, ncol=nrow(features_select), 0)
                feats <- features_select$Feature_name
                delta_T1 <- base::as.data.frame(mat)
                colnames(delta_T1) <- feats
                delta_T2 <- base::as.data.frame(mat)
                colnames(delta_T2) <- feats
                delta_M1 <- base::as.data.frame(mat)
                colnames(delta_M1) <- feats
                delta_M2 <- base::as.data.frame(mat)
                colnames(delta_M2) <- feats
            }

            else{
                graph_peaks <- list()
                peaks_quant <- list()

                # 1-num_peaks_store stores the quantification for top closest
                # peaks while last stores total quantification of the window
                for(p in 1:(num_peaks_store + 1)){
                    mat <- matrix(nrow=10, ncol=nrow(features_select), 0)
                    Intensities <- base::as.data.frame(mat)
                    colnames(Intensities) <- features_select$Feature_name

                    # If cut off pvalue is defined then store background
                    # quantifications separately
                    if(ion_intensity_cutoff == TRUE){
                        mt <- matrix(nrow=10, ncol=nrow(features_select), 0)
                        Intensities_signal_background <- base::as.data.frame(mt)
                        fn <- features_select$Feature_name
                        colnames(Intensities_signal_background) <- fn
                    }

                    # Store data for calculating scoring per sample and feature
                    mat <- matrix(nrow=1, ncol=nrow(features_select), 0)
                    feats <- features_select$Feature_name
                    delta_T1 <- base::as.data.frame(mat)
                    colnames(delta_T1) <- feats
                    delta_T2 <- base::as.data.frame(mat)
                    colnames(delta_T2) <- feats
                    delta_M1 <- base::as.data.frame(mat)
                    colnames(delta_M1) <- feats
                    delta_M2 <- base::as.data.frame(mat)
                    colnames(delta_M2) <- feats

                    isb <- Intensities_signal_background
                    peaks_quant[[p]] <- list(Intensities=Intensities,
                                             Intensities_signal_background=isb,
                                             delta_T1=delta_T1,
                                             delta_T2=delta_T2,
                                             delta_M1=delta_M1,
                                             delta_M2=delta_M2)
                }

                nms <- base::paste("Peak_", 1:num_peaks_store, sep="")
                names(peaks_quant)[1:num_peaks_store] <- nms
                names(peaks_quant)[num_peaks_store + 1] <- "Standard"
            }

            print("Extracting intensities")
            max <- nrow(features_select)
            title <- base::paste("Extracting intensities:", Sample_ID)
            label <- base::paste(round(0 / max * 100, 0), "% done")
            pb <- tcltk::tkProgressBar(title=title, label=label, min=0, max=max,
                                       width=300)
            start_time <- Sys.time()
            updatecounter <- 0
            time_require <- 0
            end <- max

            Sample_ID_save <- Sample_ID
            Sample_ID <- base::gsub("-", ".", Sample_ID)

            iso_mult <- (1:3) * 1.002054

            rtl_sel <- features_select$RT_length
            peak_width_stats <- grDevices::boxplot.stats(rtl_sel)$stats
            dif1 <- features_select$m.z - features_select$m.z_range_min
            dif2 <- features_select$RT - features_select$RT_range_min
            delta_mz <- stats::median(dif1, na.rm=TRUE)
            delta_rt <- stats::median(dif2, na.rm=TRUE) / 2

            krt <- stringr::str_split(features_select$Observed_RT, ";",
                                      simplify=TRUE)
            known_RTs <- base::as.data.frame(krt)
            loc <- which(grepl("RT_calibration", colnames(features_select)))
            colnames(known_RTs) <- base::substr(colnames(features_select)[loc],
                                                16, 1000)
            known_RTs[] <- lapply(known_RTs,
                                  function(x) as.numeric(as.character(x)))
            known_RTs_all <- known_RTs
            known_RTs <- known_RTs[, Sample_ID]

            # Expand RT window around expected window for 2Dpeak detection
            RT_window_expand_factor <- 5

            if(MassSpec_mode == "TIMSToF"){
                ikl_sel <- features_select$Inv_K0_length
                IM_width_stats <- grDevices::boxplot.stats(ikl_sel)$stats
                dif <- features_select$Inv_K0 - features_select$Inv_K0_range_min
                delta_im <- stats::median(dif, na.rm=TRUE) / 2

                oim <- stringr::str_split(features_select$Observed_IM, ";",
                                          simplify=TRUE)
                known_IMs <- base::as.data.frame(oim)
                loc <- which(grepl("IM_calibration", colnames(features_select)))
                cnames <- colnames(features_select)[loc]
                colnames(known_IMs) <- base::substr(cnames, 16, 1000)
                known_IMs[] <- lapply(known_IMs,
                                      function(x) as.numeric(as.character(x)))
                known_IMs_all <- known_IMs
                known_IMs <- known_IMs[, Sample_ID]

                IM_window_expand_factor <- 5
                dif <- features_select$Inv_K0 - features_select$Inv_K0_range_min
                delta_im <- stats::median(dif, na.rm=TRUE) / 2
            }

            if(plots == TRUE){
                dir.create(base::paste(path, "/2Dpeakselection", sep=""))
                dir.create(base::paste(path, "/2Dpeakselection/", Sample_ID,
                                       sep=""))
            }

            for(i in 1:nrow(features_select)){
                if(MassSpec_mode == "Orbitrap"){
                    # Extract all relevant ions in RT window
                    loc <- base::paste("RT_calibration.", Sample_ID, sep="")
                    RT_correction <- ifelse(RT_calibration == TRUE,
                                            features_select[i, loc], 0)

                    if(peak_detection == TRUE){
                        mn <- features_select$RT_range_min[i]
                        mx <- features_select$RT_range_max[i]
                        RT_window <- c(mn + RT_correction, mx + RT_correction)
                        ef <- (delta_rt * RT_window_expand_factor)
                        RT_window <- c(RT_window[1] - er, RT_window[2] + er)
                    }

                    else{
                        fsel_cor <- features_select$RT[i] + RT_correction
                        fsel_len <- (features_select$RT_length[i])
                        RT_window <- c(fsel_cor - fsel_len, fsel_cor + fsel_len)
                    }

                    # Get search range in spectra
                    c1 <- RT_window[1] <= max(Indices$RT_start)
                    c2 <- RT_window[2] <= max(Indices$RT_end)

                    if(c1 & c2){
                        loc1 <- which(Indices$RT_start >= RT_window[1])[1] - 1
                        loc2 <- which(Indices$RT_end >= RT_window[2])[1]
                        search_range <- (loc1):(loc2)
                    }

                    else{
                        search_range <- nrow(Indices)
                    }

                    if(any(search_range < 1)){
                        search_range <- search_range[-which(search_range < 1)]
                    }
                }

                else{
                    # Extract all relevant ions in RT window
                    loc_rt <- base::paste("RT_calibration.", Sample_ID, sep="")
                    loc_im <- base::paste("IM_calibration.", Sample_ID, sep="")
                    RT_correction <- ifelse(RT_calibration == TRUE,
                                            features_select[i, loc_rt], 0)
                    IM_correction <- features_select[i, loc_im]

                    if(peak_detection == TRUE){
                        mn <- features_select$RT_range_min[i]
                        mx <- features_select$RT_range_max[i]
                        RT_window <- c(mn + RT_correction, mx + RT_correction)
                        rtef <- (delta_rt * RT_window_expand_factor)
                        RT_window <- c(RT_window[1] - rtef, RT_window[2] + rtef)

                        mn <- features_select$Inv_K0_range_min[i]
                        mx <- features_select$Inv_K0_range_max[i]
                        IM_window <- c(mn + IM_correction, mx + IM_correction)
                        imef <- (delta_im * IM_window_expand_factor)
                        IM_window <- c(IM_window[1] - imef, IM_window[2] + imef)
                    }

                    else{
                        rtl <- features_select$RT_length[i]
                        rt_cor <- features_select$RT[i] + RT_correction
                        RT_window <- c(rt_cor - rtl, rt_cor + rtl)
                        sel_ik0 <- features_select$Inv_K0[i] + IM_correction
                        ikl <- features_select$Inv_K0_length[i]
                        IM_window <- c(sel_ik0 - ikl, sel_ik0 + ikl)
                    }

                    # Get search range in spectra
                    c1 <- RT_window[1] <= max(Indices$RT_start)
                    c2 <- RT_window[2] <= max(Indices$RT_end)

                    if(c1 & c2){
                        rw1 <- RT_window[1]
                        rw2 <- RT_window[2]
                        search_range <- which(Indices$RT_start <= rw1
                                              & Indices$RT_end >= rw1
                                              | Indices$RT_start >= rw1
                                              & Indices$RT_end <= rw2
                                              | Indices$RT_start >= rw1
                                              & Indices$RT_start <= rw2)
                        iw1 <- IM_window[1]
                        iw2 <- IM_window[2]
                        search_range_IM <- which(Indices$IM_start <= iw1
                                                 & Indices$IM_end >= iw1
                                                 | Indices$IM_start >= iw1
                                                 & Indices$IM_end <= iw2
                                                 | Indices$IM_start >= iw1
                                                 & Indices$IM_start <= iw2)
                        loc <- which(search_range %in% search_range_IM)
                        search_range <- search_range[loc]
                    }

                    else{
                        search_range <- nrow(Indices)
                    }

                    if(any(search_range < 1)){
                        search_range <- search_range[-which(search_range < 1)]
                    }
                }

                if(length(search_range) > 0){
                    sub <- data.table::rbindlist(data_frags[search_range])

                    res <- NULL
                    if(peak_detection == TRUE){
                        if(plots == TRUE){
                            fname <- features_select$Feature_name[i]
                            grDevices::pdf(base::paste(path,
                                                       "/2Dpeakselection/",
                                                       Sample_ID, "/", fname,
                                                       ".pdf", sep=""))
                        }

                        if(MassSpec_mode == "Orbitrap"){
                            res <- feature_2D_peak_selection(
                                mz_bw=0.002,
                                all_ion_data=sub,
                                selected_features=features_select[i, ],
                                cur_sample=Sample_ID,
                                known_RT=known_RTs[i],
                                delta_mz=delta_mz,
                                delta_rt=delta_rt,
                                max_delta_RT=2 * delta_rt,
                                max_delta_mz=3 * delta_mz,
                                peak_min_ion_count=peak_min_ion_count,
                                n_raster_dens_matrix=kde_resolution,
                                num_peaks_store=num_peaks_store,
                                plot=plots,
                                MassSpec_mode=MassSpec_mode
                            )
                        }

                        else{
                            res <- feature_2D_peak_selection(
                                mz_bw=0.01,
                                all_ion_data=sub,
                                selected_features=features_select[i, ],
                                cur_sample=Sample_ID,
                                known_RT=known_RTs[i],
                                delta_mz=delta_mz,
                                delta_rt=delta_rt,
                                max_delta_RT=2 * delta_rt,
                                max_delta_mz=3 * delta_mz,
                                peak_min_ion_count=peak_min_ion_count,
                                n_raster_dens_matrix=kde_resolution,
                                num_peaks_store=num_peaks_store,
                                plot=plots,
                                MassSpec_mode=MassSpec_mode,
                                delta_im=delta_im,
                                IM_window_expand_factor=IM_window_expand_factor,
                                use_IM_data=use_IM_data
                            )
                        }

                        if(plots == TRUE){
                            grDevices::dev.off()
                        }

                        res$graph <- NULL
                    }

                    # No peak detection
                    else{
                        res <- feature_no_peak_selection(
                            all_ion_data=sub,
                            selected_features=features_select[i, ],
                            cur_sample=Sample_ID,
                            num_peaks_store=num_peaks_store,
                            MassSpec_mode=MassSpec_mode,
                            use_IM_data=use_IM_data)
                    }

                    # Check if any ions are available in any window
                    temp <- unlist(res)

                    if(any(temp[grepl("ion_count", names(temp))] > 0)){
                        # No peak detection - simply sum all ion intensities
                        # within expected window
                        if(peak_detection == FALSE){
                            # All ions
                            pki <- res$Standard$Peak_info
                            m.z_window_final <- c(pki$mz_window_lower,
                                                  pki$mz_window_upper)
                            RT_window_final <- c(pki$RT_win_lower,
                                                 pki$RT_win_upper)
                            res_int_total <- res$Standard$ion_data$Intensity
                            res_mz_total <- res$Standard$ion_data$m.z
                            res_rt_total <- res$Standard$ion_data$RT
                            res_isotp_total <- res$Standard$ion_data$isotope
                            length_data_total <- length(res_int_total)

                            if(MassSpec_mode == "TIMSToF"){
                                IM_window_final <- c(pki$IM_window_lower,
                                                     pki$IM_window_upper)
                                res_im_total <- res$Standard$ion_data$`1/K0`
                            }

                            # Determine which ions are showing an intensity
                            # above the background signal intensity
                            if(ion_intensity_cutoff == TRUE){
                                ity <- res$Standard$ion_data$Intensity
                                bk_md <- mean_bkgr_ion_inty_model[i, Sample_ID]
                                sdi <- sd_bkgr_ion_inty[1, Sample_ID]
                                x <- ifelse(base::log2(ity) > bk_md + 2 * sdi,
                                            1, 0)
                                res$Standard$ion_data$signal_background <- x

                                sbk <- res$Standard$ion_data$signal_background
                                loc <- which(sbk == 1)
                                res_int_signal <- res_int_total[loc]
                                res_mz_signal <- res_mz_total[loc]
                                res_rt_signal <- res_rt_total[loc]
                                res_isotope_signal <- res_isotp_total[loc]
                                length_data_signal <- length(res_int_signal)

                                # Save signal Intensities
                                if(length(res_int_signal) > 0){
                                    logs <- log10(res_int_signal)
                                    data.table::set(
                                        Intensities_signal_background,
                                        i=as.integer(1:8),
                                        j=as.integer(i),
                                        value=list(c(log10(sum(res_int_signal)),
                                                     length(res_int_signal),
                                                     mean(logs, na.rm=TRUE),
                                                     stats::sd(logs,
                                                               na.rm=TRUE),
                                                     m.z_window_final[1],
                                                     m.z_window_final[2],
                                                     RT_window_final[1],
                                                     RT_window_final[2]))
                                    )
                                }

                                # Save signal+background Intensities
                                if(length(res_int_total) > 0){
                                    logt <- log10(res_int_total)
                                    data.table::set(
                                        Intensities_signal_background,
                                        i=as.integer(1:8),
                                        j=as.integer(i),
                                        value=list(c(log10(sum(res_int_total)),
                                                     length(res_int_total),
                                                     mean(logt, na.rm=TRUE),
                                                     stats::sd(logt,
                                                               na.rm=TRUE),
                                                     m.z_window_final[1],
                                                     m.z_window_final[2],
                                                     RT_window_final[1],
                                                     RT_window_final[2]))
                                    )
                                }
                            }

                            # No discrimination between signal and background
                            # ions. Save all ions
                            else{
                                # Save signal+background Intensities
                                if(length(res_int_total) > 0){
                                    logt <- log10(res_int_total)
                                    data.table::set(
                                        Intensities,
                                        i=as.integer(1:8),
                                        j=as.integer(i),
                                        value=list(c(log10(sum(res_int_total)),
                                                     length(res_int_total),
                                                     mean(logt, na.rm=TRUE),
                                                     stats::sd(logt,
                                                               na.rm=TRUE),
                                                     m.z_window_final[1],
                                                     m.z_window_final[2],
                                                     RT_window_final[1],
                                                     RT_window_final[2]))
                                    )
                                }
                            }

                            # Feature alignment Scoring based on algorithm from
                            # DeMix-Q (Zhang, 2016) - use all ions
                            if(length_data_total > 0){
                                df <- base::data.frame(mz=res_mz_total,
                                                       rt=res_rt_total,
                                                       iso=res_isotp_total,
                                                       int=res_int_total)
                                df_mono_isotp <- df[which(df$iso == 0), ]
                                df_isotope_1 <- df[which(df$iso == 1), ]

                                if(nrow(df_mono_isotp) > 0){
                                    # Deviation from consensus feature
                                    int <- df_mono_isotp$int
                                    max_mono <- which(int == max(int,
                                                                 na.rm=TRUE))
                                    max_rt <- df_mono_isotp$rt[max_mono]
                                    max_mz <- df_mono_isotp$mz[max_mono]
                                    frt <- (features_select$RT[i])
                                    fmz <- (features_select$m.z[i])
                                    data.table::set(delta_T1, i=1L,
                                                    j=as.integer(i),
                                                    mean(max_rt, na.rm=TRUE)
                                                    - frt)
                                    data.table::set(delta_M1, i=1L,
                                                    j=as.integer(i),
                                                    mean(max_mz, na.rm=TRUE)
                                                    - fmz)

                                    if(nrow(df_isotope_1) > 0){
                                        # Deviation from monoisotopic ion to M+1
                                        # isotope ion
                                        int <- df_isotope_1$int
                                        mint <- max(int, na.rm=TRUE)
                                        max_iso_1 <- which(int == mint)
                                        max_rt <- df_mono_isotp$rt[max_mono]
                                        max_mz <- df_mono_isotp$mz[max_mono]
                                        mrt <- mean(df_isotope_1$rt[max_iso_1],
                                                    na.rm=TRUE)
                                        mmz <-mean(df_isotope_1$mz[max_iso_1],
                                                   na.rm=TRUE)
                                        ichar <- (iso_mult[1]
                                                  / features_select$Charge[i])
                                        dif1 <- mean(max_rt, na.rm=TRUE) - mrt
                                        dif2 <- ((mean(max_mz, na.rm=TRUE)
                                                  + ichar) - mmz)
                                        data.table::set(delta_T2, i=1L,
                                                        j=as.integer(i), dif1)
                                        data.table::set(delta_M2, i=1L,
                                                        j=as.integer(i), dif2)
                                    }
                                }
                            }
                        }

                        # With peak detection
                        else{
                            #Store up to num_peaks_store potential peak windows
                            # + the original expected window
                            num_peaks <- length(which(names(res) != "Standard"
                                                      & names(res) != "graph"))
                            for(p in names(res)[which(names(res) != "graph")]){
                                pki <- res[[p]]$Peak_info
                                m.z_window_final <- c(pki$mz_window_lower,
                                                      pki$mz_window_upper)
                                RT_window_final <- c(pki$RT_win_lower,
                                                     pki$RT_win_upper)

                                res_int_total <- res[[p]]$ion_data$Intensity
                                res_mz_total <- res[[p]]$ion_data$m.z
                                res_rt_total <- res[[p]]$ion_data$RT
                                res_isotp_total <- res[[p]]$ion_data$isotope
                                length_data_total <- length(res_int_total)

                                if(MassSpec_mode == "TIMSToF"){
                                    stpki <- res$Standard$Peak_info
                                    IM_window_final <- c(stpki$IM_window_lower,
                                                         stpki$IM_window_upper)
                                    res_im_total <- res$Standard$ion_data$`1/K0`
                                }

                                # Determine which ions are showing an intensity
                                # above the background signal intensity
                                if(ion_intensity_cutoff == TRUE){
                                    ity <- res[[p]]$ion_data$Intensity
                                    logi <- base::log2(ity)
                                    mkim <- mean_bkgr_ion_inty_model[i,
                                                                     Sample_ID]
                                    sdbi <- sd_bkgr_ion_inty[1, Sample_ID]
                                    x <- ifelse(logi > mkim + 2 * sdbi, 1, 0)
                                    res[[p]]$ion_data$signal_background <- x

                                    loc <- which(x == 1)
                                    res_int_signal <- res_int_total[loc]
                                    res_mz_signal <- res_mz_total[loc]
                                    res_rt_signal <- res_rt_total[loc]
                                    res_isotope_signal <- res_isotp_total[loc]
                                    length_data_signal <- length(res_int_signal)

                                    pk <- peaks_quant[[p]]

                                    # Save signal Intensities
                                    if(length(res_int_signal) > 0){
                                        pkint <- pk$Intensities
                                        logs <- log10(res_int_signal)
                                        lst <- list(c(
                                            log10(sum(res_int_signal)),
                                            length(res_int_signal),
                                            mean(logs, na.rm=TRUE),
                                            stats::sd(logs, na.rm=TRUE),
                                            m.z_window_final[1],
                                            m.z_window_final[2],
                                            RT_window_final[1],
                                            RT_window_final[2],
                                            num_peaks,
                                            pki$known_peak
                                        ))
                                        data.table::set(pkint,
                                                        i=as.integer(1:10),
                                                        j=as.integer(i),
                                                        value=lst)
                                    }

                                    # Save signal+background Intensities
                                    if(length(res_int_total) > 0){
                                        pkib <- pk$Intensities_signal_background
                                        logi <- log10(res_int_total)
                                        lst <- list(c(log10(sum(res_int_total)),
                                                      length(res_int_total),
                                                      mean(logi, na.rm=TRUE),
                                                      stats::sd(logi,
                                                                na.rm=TRUE),
                                                      m.z_window_final[1],
                                                      m.z_window_final[2],
                                                      RT_window_final[1],
                                                      RT_window_final[2],
                                                      num_peaks,
                                                      pki$known_peak))
                                        data.table::set(pkib,
                                                        i=as.integer(1:10),
                                                        j=as.integer(i),
                                                        value=lst)
                                    }
                                    # No intensity at all then at least save
                                    # standard information
                                    else{
                                        pkib <- pk$Intensities_signal_background
                                        lst <- list(c(NA, length(res_int_total),
                                                      NA, NA,
                                                      m.z_window_final[1],
                                                      m.z_window_final[2],
                                                      RT_window_final[1],
                                                      RT_window_final[2],
                                                      num_peaks,
                                                      pki$known_peak))
                                        data.table::set(pkib,
                                                        i=as.integer(1:10),
                                                        j=as.integer(i),
                                                        value=lst)
                                    }
                                }

                                # No discrimination between signal and
                                # background ions. Save all ions
                                else{
                                    # Save signal+background Intensities
                                    if(length(res_int_total) > 0){
                                        logt <- log10(res_int_total)
                                        lst <- list(c(log10(sum(res_int_total)),
                                                      length(res_int_total),
                                                      mean(logt, na.rm=TRUE),
                                                      stats::sd(logt,
                                                                na.rm=TRUE),
                                                      m.z_window_final[1],
                                                      m.z_window_final[2],
                                                      RT_window_final[1],
                                                      RT_window_final[2],
                                                      num_peaks,
                                                      pki$known_peak))
                                        data.table::set(Intensities,
                                                        i=as.integer(1:10),
                                                        j=as.integer(i),
                                                        value=lst)
                                    }

                                    # No intensity at all then at least save
                                    # standard information
                                    else{
                                        lst <- list(c(NA, length(res_int_total),
                                                      NA, NA,
                                                      m.z_window_final[1],
                                                      m.z_window_final[2],
                                                      RT_window_final[1],
                                                      RT_window_final[2],
                                                      num_peaks,
                                                      pki$known_peak))
                                        data.table::set(Intensities,
                                                        i=as.integer(1:10),
                                                        j=as.integer(i),
                                                        value=lst)
                                    }
                                }

                                # Feature alignment Scoreing based on algorithm
                                # from DeMix-Q (Zhang, 2016) - use all ions
                                if(length_data_total > 0){
                                    df <- base::data.frame(mz=res_mz_total,
                                                           rt=res_rt_total,
                                                           iso=res_isotp_total,
                                                           int=res_int_total)
                                    df_mono_isotp <- df[which(df$iso == 0), ]
                                    df_isotope_1 <- df[which(df$iso == 1), ]

                                    if(nrow(df_mono_isotp) > 0){
                                        # Deviation from consensus feature
                                        pkq <- peaks_quant[[p]]
                                        mx <- max(df_mono_isotp$int,
                                                  na.rm=TRUE)
                                        mxmn <- which(df_mono_isotp$int == mx)
                                        mrt <- mean(df_mono_isotp$rt[mxmn],
                                                    na.rm=TRUE)
                                        frt <- (features_select$RT[i])
                                        mmz <- mean(df_mono_isotp$mz[mxmn],
                                                    na.rm=TRUE)
                                        fmz <- (features_select$m.z[i])
                                        data.table::set(pkq$delta_T1, i=1L,
                                                        j=as.integer(i),
                                                        mrt - frt)
                                        data.table::set(pkq$delta_M1, i=1L,
                                                        j=as.integer(i),
                                                        mmz - fmz)

                                        if(nrow(df_isotope_1) > 0){
                                            # Deviation from monoisotopic ion to
                                            # M+1 isotope ion
                                            mx <- max(df_isotope_1$int,
                                                      na.rm=TRUE)
                                            mx1 <- which(df_isotope_1$int == mx)
                                            chrg <- features_select$Charge[i]
                                            v1 <- (mean(df_mono_isotp$rt[mxmn],
                                                        na.rm=TRUE)
                                                   - mean(df_isotope_1$rt[mx1],
                                                          na.rm=TRUE))
                                            v2 <- ((mean(df_mono_isotp$mz[mxmn],
                                                         na.rm=TRUE)
                                                    + (iso_mult[1] / chrg))
                                                   - mean(df_isotope_1$mz[mx1],
                                                          na.rm=TRUE))
                                            data.table::set(pkq$delta_T2, i=1L,
                                                            j=as.integer(i), v1)
                                            data.table::set(pkq$delta_M2, i=1L,
                                                            j=as.integer(i), v2)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                updatecounter <- updatecounter + 1

                if(updatecounter >= 5){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d:%02d', td@day,
                                            td@hour, lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))

                    updatecounter <- 0
                    label <- base::paste(round(i / max * 100, 0), " % done (",
                                         i, "/", max, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }

            }

            close(pb)

            if(peak_detection == FALSE){
                Intensities <- t(Intensities)
                Intensities[Intensities == 0] <- NA

                if(ion_intensity_cutoff == TRUE){
                    int_sig_bkg <- t(Intensities_signal_background)
                    int_sig_bkg[int_sig_bkg == 0] <- NA
                }

                delta_T1 <- t(delta_T1)
                delta_T1[delta_T1 == 0] <- NA
                delta_M1 <- t(delta_M1)
                delta_M1[delta_M1 == 0] <- NA
                delta_T2 <- t(delta_T2)
                delta_T2[delta_T2 == 0] <- NA
                delta_M2 <- t(delta_M2)
                delta_M2[delta_M2 == 0] <- NA

                peaks_quant <- list()

                if(ion_intensity_cutoff == TRUE){
                    lst <- list(Intensities=Intensities,
                                Intensities_signal_background=int_sig_bkg,
                                delta_T1=delta_T1, delta_T2=delta_T2,
                                delta_M1=delta_M1, delta_M2=delta_M2)
                    peaks_quant[[1]] <- lst
                }

                else{
                    peaks_quant[[1]] <- list(Intensities=Intensities,
                                             delta_T1=delta_T1,
                                             delta_T2=delta_T2,
                                             delta_M1=delta_M1,
                                             delta_M2=delta_M2)
                }

                return(list(peaks_quant=peaks_quant, graph_peaks=graph_peaks))
            }

            else{
                for(p in 1:(num_peaks_store + 1)){
                    transp <- t(peaks_quant[[p]]$Intensities)
                    peaks_quant[[p]]$Intensities <- transp
                    loc <- peaks_quant[[p]]$Intensities == 0
                    peaks_quant[[p]]$Intensities[loc] <- NA

                    trp <- t(peaks_quant[[p]]$Intensities_signal_background)
                    peaks_quant[[p]]$Intensities_signal_background <- trp
                    loc <- peaks_quant[[p]]$Intensities_signal_background == 0
                    peaks_quant[[p]]$Intensities_signal_background[loc] <- NA

                    peaks_quant[[p]]$delta_T1 <- t(peaks_quant[[p]]$delta_T1)
                    loc <- peaks_quant[[p]]$delta_T1 == 0
                    peaks_quant[[p]]$delta_T1[loc] <- NA
                    peaks_quant[[p]]$delta_M1 <- t(peaks_quant[[p]]$delta_M1)
                    loc <- peaks_quant[[p]]$delta_M1 == 0
                    peaks_quant[[p]]$delta_M1[loc] <- NA
                    peaks_quant[[p]]$delta_T2 <- t(peaks_quant[[p]]$delta_T2)
                    loc <- peaks_quant[[p]]$delta_T2 == 0
                    peaks_quant[[p]]$delta_T2[loc] <- NA
                    peaks_quant[[p]]$delta_M2 <- t(peaks_quant[[p]]$delta_M2)
                    loc <- peaks_quant[[p]]$delta_M2 == 0
                    peaks_quant[[p]]$delta_M2[loc] <- NA
                }

                return(list(peaks_quant=peaks_quant, graph_peaks=graph_peaks))
            }
        }

        # Function to call get_intensities and finally save resulting table as
        # .tab data table
        extract_intensities <- function(Sample_ID, features_select, path_to_raw,
                                        path_to_output_folder, RT_calibration,
                                        mz_calibration, peak_detection,
                                        ion_intensity_cutoff,
                                        mean_bkgr_ion_inty_model,
                                        sd_bkgr_ion_inty,
                                        peak_min_ion_count, kde_resolution,
                                        num_peaks_store, plots, MassSpec_mode,
                                        use_IM_data){
            res <- get_intensities(
                Sample_ID,
                path=path_to_raw,
                features_select=features_select,
                RT_calibration=RT_calibration,
                mz_calibration=mz_calibration,
                peak_detection=peak_detection,
                ion_intensity_cutoff=ion_intensity_cutoff,
                mean_bkgr_ion_inty_model=mean_bkgr_ion_inty_model,
                sd_bkgr_ion_inty=sd_bkgr_ion_inty,
                peak_min_ion_count=peak_min_ion_count,
                kde_resolution=kde_resolution,
                num_peaks_store=num_peaks_store,
                plots=plots,
                MassSpec_mode=MassSpec_mode,
                use_IM_data=use_IM_data
            )
            peaks_quant <- res$peaks_quant

            save(peaks_quant, file=base::paste(path_to_output_folder, "/",
                                               Sample_ID,
                                               "_feature_quant.RData", sep=""))
        }

        # Prepare threads to run extraction. The task is using much memory so
        # that more than 2 threads in parallel on a pc with 16 gb of ram
        # Results in slower performance than for just 2 threads
        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)
        res <- foreach::foreach(i=Sample_IDs) %dopar% {
            extract_intensities(
                Sample_ID=i,
                features_select=features_select,
                path_to_raw,
                path_to_output_folder,
                RT_calibration,
                mz_calibration,
                peak_detection=peak_detection,
                ion_intensity_cutoff=ion_intensity_cutoff,
                mean_bkgr_ion_inty_model=mean_bkgr_ion_inty_model,
                sd_bkgr_ion_inty=sd_bkgr_ion_inty,
                peak_min_ion_count=peak_min_ion_count,
                kde_resolution=kde_resolution,
                num_peaks_store=num_peaks_store,
                plots=plots,
                MassSpec_mode=MassSpec_mode,
                use_IM_data=use_IM_data)
        }

        parallel::stopCluster(cl)
    }

    if(MassSpec_mode == "Orbitrap"){
        # Check for which samples mzXML files are available
        mzXMLfiles <- list.files(base::paste(path_to_mzXML, "/all_ion_lists",
                                             sep=""))
        mzXMLfiles <- mzXMLfiles[which(grepl(".RData", mzXMLfiles))]
        samples <- mzXMLfiles
        samples <- base::substr(samples, 1,
                                regexpr("_all_ions.RData", samples) - 1)
    }

    if(MassSpec_mode == "TIMSToF"){
        # Check for which samples extracted spectra files are available
        mzXMLfiles <- list.files(path_to_extracted_spectra)
        mzXMLfiles <- mzXMLfiles[which(grepl(".RData", mzXMLfiles))]
        samples <- mzXMLfiles
        samples <- base::substr(samples, 1,
                                regexpr("_all_ions.RData", samples) - 1)
        path_to_mzXML <- base::gsub("\\\\all_ion_lists|/all_ion_lists", "",
                                    path_to_extracted_spectra)
    }

    # Keep samples which should be actually requantified
    loc <- which(grepl("RT_calibration\\.", colnames(features)))
    Requant_samples <- base::gsub("RT_calibration\\.", "",
                                  colnames(features)[loc])
    Requant_samples <- base::gsub("\\.", "-", Requant_samples)

    samples <- samples[which(samples %in% Requant_samples)]

    # Perform quantification of decoy features
    pth <- base::paste(path_to_mzXML,
                       "/all_ion_lists/Extracted decoy intensities",
                       output_file_names_add, sep="")
    dir.create(pth, showWarnings=FALSE)
    available <- list.files(pth)
    available <- base::gsub("_feature_quant.RData", "", available)

    # Return a warning if for some samples decoy features were already extracted
    if(length(which(samples %in% available)) > 0){
        indx <- which(samples %in% available)
        options(warn=1)
        warning(paste0(Sys.time(), " Decoy intensities were already extracted ",
                       "for ", paste(samples[indx], collapse=","), " and will ",
                       "be used for subsequent IceR steps. If this is ",
                       "unintended because raw files, MaxQuant results, or ",
                       "IceR parameters were changed, please stop IceR now, ",
                       "delete the folder ", pth, " and restart IceR. ",
                       "Consider removing ", pth, " as well. If any of the ",
                       "previously mentioned changes were made, please ",
                       "consider removing the complete Temporary_files folder ",
                       "to enable a fresh run of IceR or at least ",
                       "Quantification_raw_results.RData."))
        options(warn=-1)
    }

    if(length(which(samples %not in% available)) > 0){
        print(paste0(Sys.time(), " Extract decoy intensities"))
        selected_decoys <- which(grepl("_d", features$Feature_name))

        samples <- samples[which(samples %not in% available)]
        # Single core in case of TIMSToF data as it fills memory too much
        extract_intensities_worker(
            Sample_IDs=as.character(samples),
            features_select=features[selected_decoys, ],
            path_to_raw=base::paste(path_to_mzXML, "/all_ion_lists", sep=""),
            path_to_output_folder=pth,
            RT_calibration=RT_calibration,
            mz_calibration=mz_calibration,
            peak_detection=FALSE,
            n_cores=ifelse(MassSpec_mode == "Orbitrap", n_cores,
                           ifelse(n_cores >= 3, 3, n_cores)),
            MassSpec_mode=MassSpec_mode,
            use_IM_data=use_IM_data
        )
        print(paste0(Sys.time(), " Extract decoy intensities finished"))
    }

    # Determine distribution of background ion intensities
    files <- list.files(pth)
    samples <- files[which(grepl("_feature_quant.RData", files))]
    samples <- base::substr(samples, 1,
                            regexpr("_feature_quant.RData", samples) - 1)

    if(length(which(!grepl("_feature_quant.RData", files))) > 0){
        files <- files[-which(!grepl("_feature_quant.RData", files))]
    }

    features_select <- features[which(grepl("_d", features$Feature_name)), ]

    decoy_intensities <- base::as.data.frame(matrix(ncol=length(samples),
                                                    nrow=nrow(features_select)))
    colnames(decoy_intensities) <- samples
    decoy_intensities <- sapply(decoy_intensities, as.numeric)
    decoy_intensities <- data.table::as.data.table(decoy_intensities)
    rownames(decoy_intensities) <- features_select$Feature_name

    decoy_ioncount <- decoy_intensities
    decoy_mean_intensity <- decoy_intensities
    decoy_sd_intensity <- decoy_intensities

    for(c in 1:ncol(decoy_intensities)){
        # Load stored data into variable peaks_quant
        peaks_quant <- NULL
        load(base::paste(path_to_mzXML,
                         "/all_ion_lists/Extracted decoy intensities",
                         output_file_names_add, "/",
                         colnames(decoy_intensities)[c], "_feature_quant.RData",
                         sep=""))

        signal <- peaks_quant[[1]]$Intensities

        decoy_intensities[,c] <- base::log2(10^signal[, 1])
        decoy_ioncount[, c] <- signal[, 2]
        decoy_mean_intensity[, c] <- base::log2(10^signal[, 3])
        decoy_sd_intensity[, c] <- base::log2(10^signal[, 4])
    }

    # Plot general numbers of quantifications of decoy features
    print(paste0(Sys.time(), " Fit decoy models"))
    setwd(path_to_features)
    fname <- "Temporary_files/Decoy feature quantification parameters.pdf"
    grDevices::pdf(fname)

    # Mean decoy intensity (intensity of a single decoy ion)
    RT_all <- rep(as.numeric(features_select$RT), ncol(decoy_mean_intensity))
    x_all <- as.numeric(as.matrix(decoy_mean_intensity))
    graphics::smoothScatter(RT_all, x_all, ylab="Intensity, log2",
                            main="All samples - Decoy feature mean intensity",
                            xlab="RT [min]")

    # Try to fit an average generalised additive model to determine a RT
    # dependent mean intensity and sd of intensity
    if(length(which(!is.na(unique(x_all)))) > 10){
        fit_gam_mean <- mgcv::gam(x_all ~ s(RT_all), method="REML")
    }

    else{
        warning(paste0("Too few decoy ion intensities available. Subsequent ",
                       "statistical evaluations might be wrong."))
        x_all <- rnorm(length(RT_all), mean=1, sd=0.1)
        fit_gam_mean <- mgcv::gam(x_all ~ s(RT_all), method="REML")
    }

    x_pred <- seq(min(features_select$RT, na.rm=TRUE),
                  max(features_select$RT, na.rm=TRUE),
                  length.out=nrow(features_select))
    y_pred <- stats::predict(fit_gam_mean, base::data.frame(RT_all=x_pred))
    graphics::lines(x_pred, y_pred, col="red")
    graphics::legend("topright", legend="GAM", lty=1, col="red")

    # Now fit gam models per sample
    fit_gam_per_sample <- list()

    for(c in 1:ncol(decoy_mean_intensity)){
        RT <- as.numeric(features_select$RT)
        x <- as.numeric(as.matrix(decoy_mean_intensity)[, c])
        cnames <- colnames(decoy_mean_intensity)[c]
        graphics::smoothScatter(RT, x, ylab="Mean intensity, log2",
                                main=base::paste(cnames,
                                                 "Decoy feature intensity"),
                                xlab="RT [min]")

        # Try to fit an average generalised additive model to determine a RT
        # dependent mean intensity and sd of intensity. If not enough data
        # points available, use the average gam
        gam <- tryCatch({mgcv::gam(x ~ s(RT), method="REML")},
                        error=function(error_condition){
                            RT <- rep(as.numeric(features_select$RT),
                                      ncol(decoy_mean_intensity))
                            x <- as.numeric(as.matrix(decoy_mean_intensity))
                            mgcv::gam(x ~ s(RT), method="REML")
                        })
        x_pred <- seq(min(features_select$RT, na.rm=TRUE),
                      max(features_select$RT, na.rm=TRUE),
                      length.out=nrow(features_select))
        y_pred <- stats::predict(gam, base::data.frame(RT=x_pred))
        graphics::lines(x_pred, y_pred, col="red")
        graphics::legend("topright", legend="GAM", lty=1, col="red")
        fit_gam_per_sample[[colnames(decoy_mean_intensity)[c]]] <- gam
    }

    graphics::par(mfrow=c(2, 2))
    graphics::boxplot(as.numeric(as.matrix(decoy_intensities)), outline=FALSE,
                      ylab="Summed intensity, log2",
                      main="Summed intensity of decoy ions")
    graphics::boxplot(as.numeric(as.matrix(decoy_mean_intensity)),
                      outline=FALSE, ylab="Mean intensity, log2",
                      main="Mean intensity of decoy ions")
    graphics::boxplot(as.numeric(as.matrix(decoy_sd_intensity)), outline=FALSE,
                      ylab="SD of intensity, log2",
                      main="SD of intensity of decoy ions")
    graphics::boxplot(as.numeric(as.matrix(decoy_ioncount)), outline=FALSE,
                      ylab="Count",
                      main="Number of ions in decoys with quantification")
    graphics::par(mfrow=c(1, 1))
    grDevices::dev.off()

    # Now predict background intensities per sample and feature by using
    # individually fitted GAMs, multiply with median decoy ion count and add
    # some noise using observed decoy intensity sds per sample

    # Used for determining which ions are background ions and which are signal
    # ions
    backgr_inty_GAM_tab_x_feat <- NULL
    # Used to later impute missing values
    backgr_inty_GAM_tab_x_feat_sum <- NULL

    med_sd <- matrixStats::colMedians(as.matrix(decoy_sd_intensity), na.rm=TRUE)
    median_decoy_sd_per_sample <- base::as.data.frame(t(med_sd))
    colnames(median_decoy_sd_per_sample) <- colnames(decoy_sd_intensity)

    med_ion <- matrixStats::colMedians(as.matrix(decoy_ioncount), na.rm=TRUE)
    median_decoy_ion_count_per_sample <- base::as.data.frame(t(med_ion))
    colnames(median_decoy_ion_count_per_sample) <- colnames(decoy_ioncount)

    max <- ncol(decoy_mean_intensity)
    pb <- tcltk::tkProgressBar(title="Model background intensities per feature",
                               label=base::paste(round(0 / max * 100, 0),
                               "% done"), min=0, max=max, width=300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    for(c in 1:ncol(decoy_mean_intensity)){
        cur_sample <- colnames(decoy_mean_intensity)[c]
        pred <- stats::predict(fit_gam_per_sample[[cur_sample]],
                               base::data.frame(RT=features$RT))
        df <- base::data.frame(mean_intensity=pred)

        if(c != 1){
            bndt <- cbind(backgr_inty_GAM_tab_x_feat, df)
            backgr_inty_GAM_tab_x_feat <- bndt
            bnds <- cbind(backgr_inty_GAM_tab_x_feat_sum, df)
            backgr_inty_GAM_tab_x_feat_sum <- bnds
        }

        else{
            backgr_inty_GAM_tab_x_feat <- df
            backgr_inty_GAM_tab_x_feat_sum <- df
        }

        nct <- ncol(backgr_inty_GAM_tab_x_feat)
        ncs <- ncol(backgr_inty_GAM_tab_x_feat_sum)
        colnames(backgr_inty_GAM_tab_x_feat)[nct] <- cur_sample
        colnames(backgr_inty_GAM_tab_x_feat_sum)[ncs] <- cur_sample

        rnm <- stats::rnorm(n=nrow(backgr_inty_GAM_tab_x_feat_sum),
                            mean=backgr_inty_GAM_tab_x_feat_sum[, c],
                            sd=median_decoy_sd_per_sample[1, c])
        backgr_inty_GAM_tab_x_feat_sum[, c] <- rnm
        updatecounter <- updatecounter + 1

        if(updatecounter >= 1){
            time_elapsed <- difftime(Sys.time(), start_time, units="secs")
            time_require <- (time_elapsed / (c / max)) * (1 - (c / max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                    lubridate::minute(td),
                                    round(lubridate::second(td), digits=0))
            updatecounter <- 0
            label <- base::paste(round(c / max * 100, 0), " % done (", c, "/",
                                 max, ", Time require: ", time_require, ")",
                                 sep="")
            tcltk::setTkProgressBar(pb, c, label=label)
        }
    }

    close(pb)

    sd_background_intensity <- median_decoy_sd_per_sample
    rownames(backgr_inty_GAM_tab_x_feat) <- features$Feature_name
    rownames(backgr_inty_GAM_tab_x_feat_sum) <- features$Feature_name

    df <- base::data.frame(RT_all=RT_all, mean_intensity=x_all)
    bigtxf <- backgr_inty_GAM_tab_x_feat
    bigtxfs <- backgr_inty_GAM_tab_x_feat_sum
    lst <- list(
        Decoy_mean_intensity_per_RT=df,
        GAM_model_average=fit_gam_mean,
        GAM_model_per_sample=fit_gam_per_sample,
        backgr_inty_GAM_tab_x_feat=bigtxf,
        backgr_inty_GAM_tab_x_feat_sum=bigtxfs,
        sd_background_intensity=sd_background_intensity,
        median_decoy_ion_count_per_sample=median_decoy_ion_count_per_sample,
        decoy_intensities=decoy_intensities,
        decoy_mean_intensity=decoy_mean_intensity,
        decoy_sd_intensity=decoy_sd_intensity,
        decoy_ioncount=decoy_ioncount
    )
    QC_data[["Decoy_feature_parameters"]] <- lst

    # Now select 1000 decoy features randomly for which we will perform peak
    # selection and quantification
    decoy_features <- which(grepl("_d", features$Feature_name))
    if(length(decoy_features) > 1000){
        set.seed(1)
        select_decoys <- sample(decoy_features, 1000)
        decoy_features_keep <- features[select_decoys, ]
        # Now add isotope peaks of these decoy features
        isotope_features <- decoy_features_keep
        dif <- (isotope_features$m.z * isotope_features$Charge)
        isotope_features$m.z <- (dif + 1.002054) / isotope_features$Charge
        delta_mz <- (isotope_features$m.z_range_max
                     - isotope_features$m.z_range_min)
        isotope_features$m.z_range_max <- isotope_features$m.z + (delta_mz / 2)
        isotope_features$m.z_range_min <- isotope_features$m.z - (delta_mz / 2)
        isotope_features$Feature_name <- base::paste(
            isotope_features$Feature_name,
            "_i",
            sep=""
        )
        decoy_features_keep <- rbind(decoy_features_keep, isotope_features)
    }

    # Less than 1000 decoys so just continue
    else{
        decoy_features_keep <- features[decoy_features, ]
        # Now add isotope peaks of these decoy features
        isotope_features <- decoy_features_keep
        dif <- (isotope_features$m.z * isotope_features$Charge)
        isotope_features$m.z <- (dif + 1.002054) / isotope_features$Charge
        delta_mz <- (isotope_features$m.z_range_max
                     - isotope_features$m.z_range_min)
        isotope_features$m.z_range_max <- isotope_features$m.z + (delta_mz / 2)
        isotope_features$m.z_range_min <- isotope_features$m.z - (delta_mz / 2)
        isotope_features$Feature_name <- base::paste(
            isotope_features$Feature_name,
            "_i",
            sep=""
        )
        decoy_features_keep <- rbind(decoy_features_keep, isotope_features)
    }

    # Remove all decoy features
    features <- features[-decoy_features, ]
    # Add selected decoy features plus isotope decoy features
    features <- rbind(features, decoy_features_keep)

    # See if already samples were converted
    samples <- mzXMLfiles
    samples <- base::substr(samples, 1, regexpr("_all_ions.RData", samples) - 1)
    loc <- which(grepl("RT_calibration\\.", colnames(features)))
    Requant_samples <- base::gsub("RT_calibration\\.", "",
                                  colnames(features)[loc])
    Requant_samples <- base::gsub("\\.", "-", Requant_samples)

    samples <- samples[which(samples %in% Requant_samples)]
    pth <- base::paste(path_to_mzXML,
                       "/all_ion_lists/Extracted feature intensities",
                       output_file_names_add, sep="")
    dir.create(pth, showWarnings=FALSE)
    available <- list.files(pth)
    available <- available[which(grepl("feature_quant.RData", available))]
    available <- base::gsub("_feature_quant.RData", "", available)
    dec_ioc <- as.numeric(as.matrix(decoy_ioncount))
    peak_min_ion_count <- grDevices::boxplot.stats(dec_ioc)$stats[4]

    # Return a warning if for some samples decoy features were already extracted
    if(length(which(samples %in% available)) > 0){
        indx <- which(samples %in% available)
        options(warn=1)
        warning(paste0(Sys.time(), " Feature intensities were already ",
                       "extracted for ", paste(samples[indx], collapse=","),
                       " and will be used for subsequent IceR steps. If this ",
                       "is unintended because raw files, MaxQuant results, or ",
                       "IceR parameters were changed, please stop IceR now, ",
                       "delete the folder ", pth, " and restart IceR. If ",
                       "any of the previously mentioned changes were made, ",
                       "please consider removing the complete Temporary_files ",
                       "folder to enable a fresh run of IceR or at least ",
                       "Quantification_raw_results.RData."))
        options(warn=-1)
    }

    if(length(which(samples %not in% available)) > 0){
        print(paste0(Sys.time(), " Perform peak detection"))
        samples <- samples[which(samples %not in% available)]
        # Use decoy defined cut of to distinguish background intensity from
        # signal intensity per feature
        extract_intensities_worker(
            Sample_IDs=as.character(samples),
            features_select=features,
            path_to_raw=base::paste(path_to_mzXML, "/all_ion_lists", sep=""),
            path_to_output_folder=pth,
            RT_calibration=RT_calibration,
            mz_calibration=mz_calibration,
            peak_detection=peak_detection,
            n_cores=ifelse(MassSpec_mode == "Orbitrap", n_cores,
                           ifelse(n_cores >= 3, 3, n_cores)),
            ion_intensity_cutoff=TRUE,
            mean_bkgr_ion_inty_model=backgr_inty_GAM_tab_x_feat,
            sd_bkgr_ion_inty=sd_background_intensity,
            peak_min_ion_count=peak_min_ion_count,
            kde_resolution=kde_resolution,
            num_peaks_store=num_peaks_store,
            plots=plot_peak_detection,
            MassSpec_mode=MassSpec_mode,
            use_IM_data=use_IM_data
        )
        # Define background peaks during 2DKDE as peaks with <= 75% quantile
        print(paste0(Sys.time(), " Peak detection finished"))
    }

    # Summarize data for all features and samples
    # Available sample data
    files <- list.files(pth)
    samples <- files[which(grepl("feature_quant.RData", files))]
    samples <- base::substr(samples, 1,
                            regexpr("_feature_quant.RData", samples) - 1)

    if(length(which(!grepl("feature_quant.RData", files))) > 0){
        files <- files[-which(!grepl("feature_quant.RData", files))]
    }

    features_select <- features
    mat <- matrix(ncol=length(samples), nrow=nrow(features_select))
    features_intensity <- base::as.data.frame(mat)
    colnames(features_intensity) <- samples
    features_intensity <- sapply(features_intensity, as.numeric)
    features_intensity <- data.table::as.data.table(features_intensity)
    rownames(features_intensity) <- features_select$Feature_name
    Icount_feat_sample_mat <- features_intensity

    feat_w_backgr_int <- features_intensity
    Icount_feat_w_bkgr_int <- features_intensity

    # Prepare matrices to store data
    if(peak_detection == TRUE){
        peak_quant <- list()

        for(p in 1:(num_peaks_store + 1)){
            feat_int_tmp <- features_intensity
            feat_w_backgr_int_tmp <- feat_w_backgr_int
            Ioc_feat_smp_mat_tmp <- Icount_feat_sample_mat
            Ioc_ft_w_backgr_int_tmp <- Icount_feat_w_bkgr_int
            Peak_rt_temp <- features_intensity
            Peak_mz_temp <- features_intensity
            num_peaks_temp <- features_intensity
            correct_peak_temp <- features_intensity
            Peak_rt_with_background_temp <- features_intensity
            Peak_mz_with_background_temp <- features_intensity
            num_peaks_with_background_temp <- features_intensity
            correct_peak_w_background_temp <- features_intensity

            peak_quant[[p]] <- list(
                features_intensity=feat_int_tmp,
                feat_w_backgr_int=feat_w_backgr_int_tmp,
                Icount_feat_sample_mat=Ioc_feat_smp_mat_tmp,
                Icount_feat_w_bkgr_int=Ioc_ft_w_backgr_int_tmp,
                Peak_rt=Peak_rt_temp,
                Peak_mz=Peak_mz_temp,
                num_peaks=num_peaks_temp,
                correct_peak=correct_peak_temp,
                Peak_rt_with_background=Peak_rt_with_background_temp,
                Peak_mz_with_background=Peak_mz_with_background_temp,
                num_peaks_with_background=num_peaks_with_background_temp,
                correct_peak_w_background=correct_peak_w_background_temp
            )
        }

        names(peak_quant) <- c(base::paste("Peak_", 1:(num_peaks_store),
                                           sep=""), "Standard")
    }

    max <- ncol(features_intensity)
    label <- base::paste(round(0 / max * 100, 0), "% done")
    pb <- tcltk::tkProgressBar(title="Merge quantification results",
                               label=label, min=0, max=max, width=300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    # Read the data
    for(c in 1:ncol(features_intensity)){
        if(peak_detection == FALSE){
            # Now extract intensities:
            # row1 = Intensity summed
            # row2 = num ions
            # row3 = mean intensity
            # row4 = sd intensity
            # row5 = detected optimum mz window minimum
            # row6 = detected optimum mz window maximum
            # row7 = detected optimum RT window minimum
            # row8 = detected optimum RT window maximum
            # row9 = number of detected peaks
            # row10 = index which of the peaks should be the correct one (based
            #         on known observed RT for the respective feature)

            # Load stored data into variable peaks_quant
            peaks_quant <- NULL
            load(base::paste(path_to_mzXML,
                             "/all_ion_lists/Extracted feature intensities",
                             output_file_names_add, "/",
                             colnames(features_intensity)[c],
                             "_feature_quant.RData", sep=""))

            signal <- peaks_quant[[1]]$Intensities
            signal_background <- peaks_quant[[1]]$Intensities_signal_background

            features_intensity[, c] <- signal[, 1]
            Icount_feat_sample_mat[, c] <- signal[, 2]

            feat_w_backgr_int[, c] <- signal_background[, 1]
            Icount_feat_w_bkgr_int[, c] <- signal_background[, 2]
        }

        else{
            # Load stored data into variable peaks_quant
            load(base::paste(path_to_mzXML,
                             "/all_ion_lists/Extracted feature intensities",
                             output_file_names_add, "/",
                             colnames(features_intensity)[c],
                             "_feature_quant.RData", sep=""))

            for(p in 1:(num_peaks_store + 1)){
                # Now extract intensities
                # row1 = Intensity summed
                # row2 = num ions
                # row3 = mean intensity
                # row4 = sd intensity
                # row5 = detected optimum mz window minimum
                # row6 = detected optimum mz window maximum
                # row7 = detected optimum RT window minimum
                # row8 = detected optimum RT window maximum
                # row9 = number of detected peaks
                # row10 = index which of the peaks should be the correct one
                #         (based on known observed RT for the respective
                #         feature)

                signal <- peaks_quant[[p]]$Intensities
                bk <- peaks_quant[[p]]$Intensities_signal_background

                peak_quant[[p]]$features_intensity[, c] <- signal[, 1]
                peak_quant[[p]]$Icount_feat_sample_mat[, c] <- signal[, 2]
                peak_quant[[p]]$Peak_rt[, c] <- (signal[, 7] + signal[, 8]) / 2
                peak_quant[[p]]$Peak_mz[, c] <- (signal[, 5] + signal[, 6]) / 2
                peak_quant[[p]]$num_peaks[, c] <- signal[, 9]
                peak_quant[[p]]$correct_peak[, c] <- signal[, 10]

                rt <- (bk[, 7] + bk[, 8]) / 2
                mz <- (bk[, 5] + bk[, 6]) / 2
                peak_quant[[p]]$feat_w_backgr_int[, c] <- bk[, 1]
                peak_quant[[p]]$Icount_feat_w_bkgr_int[, c] <- bk[, 2]
                peak_quant[[p]]$Peak_rt_with_background[, c] <- rt
                peak_quant[[p]]$Peak_mz_with_background[, c] <- mz
                peak_quant[[p]]$num_peaks_with_background[, c] <- bk[, 9]
                peak_quant[[p]]$correct_peak_w_background[, c] <- bk[, 10]
            }
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 1){
            time_elapsed <- difftime(Sys.time(), start_time,units="secs")
            time_require <- (time_elapsed / (c / max)) * (1 - (c / max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                    lubridate::minute(td),
                                    round(lubridate::second(td), digits=0))
            updatecounter <- 0
            tcltk::setTkProgressBar(pb, c,
                                    label=base::paste(round(c / max * 100, 0),
                                    " % done (", c, "/", max,
                                    ", Time require: ", time_require, ")",
                                    sep=""))
        }
    }

    close(pb)

    # If peak detection was performed, then check if peak selection was already
    # done. If yes, load stored data, if not, perform peak selection
    pth <- base::paste(path_to_features,
                       "/Temporary_files/Quantification_raw_results.RData",
                       sep="")

    if(file.exists(pth)){
        options(warn=1)
        warning(paste0("Peak selection was already performed and stored ",
                       "results will be used. If this is unintended because ",
                       "changed, please remove Quantification_raw_results",
                       ".RData in Raw files, MaxQuant results, or IceR ",
                       "parameters were Temporary_files folder or the ",
                       "complete Temporary_files folder to enable a fresh run ",
                       "of IceR."))
        options(warn=-1)

        load(pth)
        crap <- gc(FALSE)
    }

    else{
        print(paste0(Sys.time(), " Perform peak selection"))
        # If no peak detection was done, directly continue
        # Otherwise select per feature best peak
        if(peak_detection == TRUE){
            features_intensity <- as.matrix(features_intensity)
            Icount_feat_sample_mat <- as.matrix(Icount_feat_sample_mat)
            feat_w_backgr_int <- as.matrix(feat_w_backgr_int)
            Icount_feat_w_bkgr_int <- as.matrix(Icount_feat_w_bkgr_int)

            # Prepare matrix in which we store decisions which peaks where used
            # for quantification
            peak_selected <- matrix(ncol=ncol(feat_w_backgr_int),
                                    nrow=nrow(feat_w_backgr_int))
            colnames(peak_selected) <- colnames(feat_w_backgr_int)

            # Transfer all quantification data from Standard for quantifications
            # where no peak was detected
            mat1 <- as.matrix(peak_quant$Standard$num_peaks_with_background)
            mat2 <- as.matrix(peak_quant$Peak_1$num_peaks_with_background)
            sel <- which(is.na(mat1) & is.na(mat2))
            fint <- as.matrix(peak_quant$Standard$features_intensity)
            features_intensity[sel] <- fint[sel]
            icount <- as.matrix(peak_quant$Standard$Icount_feat_sample_mat)
            Icount_feat_sample_mat[sel] <- icount[sel]
            bkint <- as.matrix(peak_quant$Standard$feat_w_backgr_int)
            feat_w_backgr_int[sel] <- bkint[sel]
            icountbk <- as.matrix(peak_quant$Standard$Icount_feat_w_bkgr_int)
            Icount_feat_w_bkgr_int[sel] <- icountbk[sel]
            peak_selected[sel] <- num_peaks_store + 1 # Standard window

            # Next, transfer all quantification data from Peak1 for
            # quantifications where correct peak is known
            pk <- peak_quant$Peak_1
            sel <- which(as.matrix(pk$correct_peak_w_background == 1))
            features_intensity[sel] <- as.matrix(pk$features_intensity)[sel]
            aux <- as.matrix(pk$Icount_feat_sample_mat)
            Icount_feat_sample_mat[sel] <- aux[sel]
            feat_w_backgr_int[sel] <- as.matrix(pk$feat_w_backgr_int)[sel]
            aux <- as.matrix(pk$Icount_feat_w_bkgr_int)
            Icount_feat_w_bkgr_int[sel] <- aux[sel]
            peak_selected[sel] <- 1 # Closest peak = peak 1

            # If no peak is know, change from NA to 0
            loc <- is.na(peak_quant$Peak_1$correct_peak_w_background)
            peak_quant$Peak_1$correct_peak[loc] <- 0
            peak_quant$Peak_1$correct_peak_w_background[loc] <- 0

            features_intensity <- base::as.data.frame(features_intensity)
            Icount_feat_sample_mat <- base::as.data.frame(
                Icount_feat_sample_mat
            )
            feat_w_backgr_int <- base::as.data.frame(feat_w_backgr_int)
            Icount_feat_w_bkgr_int <- base::as.data.frame(
                Icount_feat_w_bkgr_int
            )
            peak_selected <- base::as.data.frame(peak_selected)

            # For the next step of peak quantifications where true peak is not
            # known we will need RT and mz correction factors per feature
            # Extract RT and mz correction factors per sample and feature
            indices_RT_correction <- which(grepl("RT_calibration",
                                                 colnames(features_select)))
            indices_mz_correction <- which(grepl("mz_calibration",
                                                 colnames(features_select)))
            cnames <- colnames(features_select)[indices_RT_correction]
            ordering_indices <- match(base::gsub("-", ".", samples),
                                      base::gsub("RT_calibration\\.", "",
                                                 cnames))
            indices_RT_correction <- indices_RT_correction[ordering_indices]
            indices_mz_correction <- indices_mz_correction[ordering_indices]

            RT_correction_factors <- features_select[, indices_RT_correction]
            mz_correction_factors <- features_select[, indices_mz_correction]
            colnames(RT_correction_factors) <- samples
            rownames(RT_correction_factors) <- features_select$Feature_name
            colnames(mz_correction_factors) <- samples
            rownames(mz_correction_factors) <- features_select$Feature_name

            # General RT and mz windows
            fsel <- features_select
            delta_mz <- stats::median(fsel$m.z - fsel$m.z_range_min, na.rm=TRUE)
            delta_rt <- stats::median(fsel$RT - fsel$RT_range_min,
                                      na.rm=TRUE) / 2

            # Split up quant results into chunks of 50000 features
            # Perform individually peak selection per chunk
            # Combine results afterwards again
            # Purpose: keep peak quant file size low otherwise memory can be
            # quickly overloaded in case of a large experiment
            chunk_size <- 50000
            num_chunks <- ceiling(nrow(features) / chunk_size)

            for(chunk in 1:num_chunks){
                start <- (chunk - 1) * chunk_size + 1
                end <- (chunk) * chunk_size

                if(end > nrow(features)){
                    end <- nrow(features)
                }

                features_temp <- features_select[start:end, ]
                peak_quant_temp <- peak_quant

                for(p in 1:(num_peaks_store + 1)){
                    for(p2 in 1:length(peak_quant_temp[[p]])){
                        pk <- peak_quant_temp[[p]][[p2]]
                        peak_quant_temp[[p]][[p2]] <- pk[start:end, ]
                    }
                }

                # Perform peak decision
                cl <- parallel::makeCluster(n_cores)
                doParallel::registerDoParallel(cl)

                res <- foreach::foreach(s=1:length(samples)) %dopar% {
                    fint <- features_intensity[start:end, s ,drop=FALSE]
                    cntfs <- Icount_feat_sample_mat[start:end, s, drop=FALSE]
                    fbki <- feat_w_backgr_int[start:end, s, drop=FALSE]
                    cfbi <- Icount_feat_w_bkgr_int[start:end, s, drop=FALSE]
                    pks <- peak_selected[start:end, s ,drop=FALSE]
                    rtcf <- RT_correction_factors[start:end, ]
                    mzcf <- mz_correction_factors[start:end, ]
                    peak_decision(
                        features_select=features_temp,
                        peak_quant=peak_quant_temp,
                        samples,
                        s=s,
                        RT_correction_factors=rtcf,
                        mz_correction_factors=mzcf,
                        features_intensity_sample=fint,
                        Ioncount_sample=cntfs,
                        feat_w_bkgrnd_inty_samp=fbki,
                        Icount_w_bkgr_smp=cfbi,
                        peak_selected_sample=pks,
                        delta_mz=delta_mz,
                        delta_rt=delta_rt,
                        peak_min_ion_count=peak_min_ion_count,
                        chunk=chunk,
                        num_chunks=num_chunks
                    )
                }
                parallel::stopCluster(cl)

                # Merge results from all threads
                for(c in 1:length(samples)){
                    data.table::set(features_intensity, as.integer(start:end),
                                    c, res[[c]]$features_intensity_sample)
                    data.table::set(Icount_feat_sample_mat,
                                    as.integer(start:end), c,
                                    res[[c]]$Ioncount_sample)
                    data.table::set(feat_w_backgr_int, as.integer(start:end), c,
                                    res[[c]]$feat_w_bkgrnd_inty_samp)
                    data.table::set(Icount_feat_w_bkgr_int,
                                    as.integer(start:end), c,
                                    res[[c]]$Icount_w_bkgr_smp)
                    data.table::set(peak_selected, as.integer(start:end), c,
                                    res[[c]]$peak_selected_sample)
                }
            }

            rownames(peak_selected) <- features$Feature_name
        }

        else{
            peak_selected <- NULL
        }

        print(paste0(Sys.time(), " Peak selection finished"))

        rownames(features_intensity) <- features$Feature_name
        rownames(Icount_feat_sample_mat) <- features$Feature_name
        rownames(feat_w_backgr_int) <- features$Feature_name
        rownames(Icount_feat_w_bkgr_int) <- features$Feature_name

        setwd(path_to_features)
        dir.create("Temporary_files")
        setwd(base::paste(path_to_features, "/Temporary_files", sep=""))
        print(paste0(Sys.time(), " Save peak detection and selection results"))
        save(features, features_intensity, Icount_feat_sample_mat,
             feat_w_backgr_int, Icount_feat_w_bkgr_int, peak_selected,
             peak_quant, QC_data, path_to_features, path_to_MaxQ_output,
             samples, delta_rt, delta_mz, peak_min_ion_count, num_peaks_store,
             RT_calibration, mz_calibration, abundance_estimation_correction,
             Quant_pVal_cut, align_var_score_thr, align_score_thr,
             mono_iso_alignment_cutoff, calc_peptide_LFQ, calc_protein_LFQ,
             MassSpec_mode, file="Quantification_raw_results.RData")

        crap <- gc(FALSE)
    }

    # Read previously generated outputs
    print(paste0(Sys.time(),
                 " Perform scoring of alignment and quantification"))
    setwd(path_to_features)
    grDevices::pdf("Temporary_files/Alignment and quantification scores.pdf")

    features$target_decoy <- ifelse(grepl("_d", features$Feature_name), "decoy",
                                    "target")

    # Tag decoy features which overlap with real features
    len <- length(which(features$target_decoy == "decoy"))
    remove_decoy_outlier <- vector(mode="logical", length=len)
    features_target <- features[which(features$target_decoy == "target"),]
    decoy_indices <- which(features$target_decoy == "decoy")

    max <- length(decoy_indices)
    pb <- tcltk::tkProgressBar(title="Detect target-decoy overlapping features",
                               label=base::paste(round(0 / max * 100, 0),
                                                 "% done"), min=0, max=max,
                               width=300)
    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    for(i in 1:length(decoy_indices)){
        index <- decoy_indices[i]
        c1 <- features_target$m.z_range_min > features$m.z_range_min[index]
        c2 <- features_target$m.z_range_min < features$m.z_range_max[index]
        c3 <- features_target$m.z_range_max > features$m.z_range_min[index]
        c4 <- features_target$m.z_range_max < features$m.z_range_max[index]
        overlap <- which(c1 & c2 | c3 & c4)

        c1 <- (features_target$RT_range_min[overlap]
               > features$RT_range_min[index])
        c2 <- (features_target$RT_range_min[overlap]
               < features$RT_range_max[index])
        c3 <- (features_target$RT_range_max[overlap]
               > features$RT_range_min[index])
        c4 <- (features_target$RT_range_max[overlap]
               < features$RT_range_max[index])
        overlap <- which(c1 & c2 | c3 & c4)

        if(length(overlap) > 0){
            remove_decoy_outlier[i] <- TRUE
        }

        updatecounter <- updatecounter + 1

        if(updatecounter >= 10){
            time_elapsed <- difftime(Sys.time(), start_time,units="secs")
            time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                    lubridate::minute(td),
                                    round(lubridate::second(td), digits=0))
            updatecounter <- 0
            tcltk::setTkProgressBar(pb, i,
                                    label=base::paste(round(i / max * 100, 0),
                                    " % done (", i, "/", max,
                                    ", Time require: ", time_require, ")",
                                    sep=""))
        }
    }

    close(pb)
    features$overlapping_decoys <- FALSE
    features$overlapping_decoys[decoy_indices] <- remove_decoy_outlier

    # t-test based determination of background quantifications
    # Determine if ion count significantly deviates from decoy feature ion
    # counts. Perform this for signal intensities and background+signal
    # intensities. Signal+Background quantification
    loc <- which(features$target_decoy == "decoy"
                 & !grepl("_d_i", features$Feature_name)
                 & features$overlapping_decoys == FALSE)
    decoy_ion_count <- as.numeric(as.matrix(Icount_feat_w_bkgr_int[loc, ]))
    decoy_ion_count[is.na(decoy_ion_count)] <- 0

    # Add a count of 1 for log scale
    decoy_ion_count <- decoy_ion_count + 1

    # Determine mean and sd of log2 transformed decoy ion counts for z scoring
    mean_decoy_count <- mean(base::log2(decoy_ion_count))
    sd_decoy_count <- stats::sd(base::log2(decoy_ion_count))

    # Prepare for testing if observed ion counts of target features is
    # significantly higher than for decoy ions (just random peak selection)
    target_ion_counts <- Icount_feat_w_bkgr_int
    rows <- nrow(target_ion_counts)
    rownames(target_ion_counts) <- rownames(Icount_feat_w_bkgr_int)
    target_ion_counts <- as.matrix(target_ion_counts)
    # Add 1 for log scale
    target_ion_counts <- target_ion_counts + 1
    zval <- (base::log2(target_ion_counts) - mean_decoy_count) / sd_decoy_count
    zscores <- base::as.data.frame(zval)

    zscore_to_pval <- function(z, alternative="greater"){
        if(alternative == "greater"){
            pval <- stats::pnorm(-z)
        }

        else if(alternative == "two.sided"){
            pval <- stats::pnorm(-abs(z))
            pval <- 2 * pval
        }

        else if(alternative == "less"){
            pval <- stats::pnorm(z)
        }

        return(pval)
    }

    # Perform significance test
    pval_sig_w_bckgrnd_quant <- base::as.data.frame(apply(zscores, 1:2,
                                                          zscore_to_pval))
    temp_pval_quant <- pval_sig_w_bckgrnd_quant
    temp_pval_quant[is.na(temp_pval_quant)] <- 1

    QC_data[["Decoy_ion_counts"]] <- decoy_ion_count
    QC_data[["Target_ion_counts"]] <- target_ion_counts
    QC_data[["Quant_pval"]] <- pval_sig_w_bckgrnd_quant

    # Plot pvalue of quantification per sample for target features
    # Only 15 samples can be plotted in one plot
    if(ncol(temp_pval_quant) <= 15){
        graphics::par(mar=c(10, 4, 4, 4))
        loc <- which(!grepl("_d", rownames(temp_pval_quant)))
        graphics::boxplot(-log10(temp_pval_quant[loc, ]), outline=FALSE,
                          main="Quantification pVals per target feature",
                          ylab="pValue, -log10",
                          names=colnames(temp_pval_quant), las=2)
        graphics::abline(h=-log10(Quant_pVal_cut), lty=2, col="red")
    }

    else{
        graphics::par(mar=c(10, 4, 4, 4))
        pages <- ceiling(ncol(temp_pval_quant) / 15)

        for(p in 1:pages){
            columns <- (((p - 1) * 15) + 1):(p * 15)
            columns <- columns[which(columns <= ncol(temp_pval_quant))]
            loc <- which(!grepl("_d", rownames(temp_pval_quant)))
            graphics::boxplot(-log10(temp_pval_quant[loc, columns]),
                              outline=FALSE,
                              main="Quantification pVals per target feature",
                              ylab="pValue, -log10",
                              names=colnames(temp_pval_quant)[columns], las=2)
            graphics::abline(h=-log10(Quant_pVal_cut), lty=2, col="red")
        }
    }

    # Plot pvalue of quantification per sample for decoy features
    # Only 15 samples can be plotted in one plot
    if(ncol(temp_pval_quant) <= 15){
        graphics::par(mar=c(10, 4, 4, 4))
        loc <- which(grepl("_d$", rownames(temp_pval_quant)))
        graphics::boxplot(-log10(temp_pval_quant[loc, ]), outline=FALSE,
                          main="Quantification pVals per decoy feature",
                          ylab="pValue, -log10",
                          names=colnames(temp_pval_quant), las=2)
        graphics::abline(h=-log10(Quant_pVal_cut), lty=2, col="red")
    }

    else{
        graphics::par(mar=c(10, 4, 4, 4))
        pages <- ceiling(ncol(temp_pval_quant) / 15)

        for(p in 1:pages){
            columns <- (((p - 1) * 15) + 1):(p * 15)
            columns <- columns[which(columns <= ncol(temp_pval_quant))]
            loc <- which(grepl("_d$", rownames(temp_pval_quant)))
            graphics::boxplot(-log10(temp_pval_quant[loc, columns]),
                              outline=FALSE,
                              main="Quantification pVals per decoy feature",
                              ylab="pValue, -log10",
                              names=colnames(temp_pval_quant)[columns], las=2)
            graphics::abline(h=-log10(Quant_pVal_cut), lty=2, col="red")
        }
    }

    # Plot per sample how many target feature quantifications were significant
    # and not significant
    loc <- which(!grepl("_d|_i", rownames(temp_pval_quant)))
    plot_data <- rbind(colSums(temp_pval_quant[loc, ] < Quant_pVal_cut,
                               na.rm=TRUE),
                       colSums(temp_pval_quant[loc,] > Quant_pVal_cut,
                               na.rm=TRUE))
    rownames(plot_data) <- c("certain", "uncertain")
    title <- "Certainty of ion accumulation for target features"

    # Only 15 samples can be plotted in one plot
    if(ncol(plot_data) <= 15){
        p <- Barplotsstacked(plot_data, AvgLine=FALSE, col=c("chocolate",
                                                             "grey"),
                             shownumbers_total=FALSE, shownumbers=TRUE,
                             Legends=rownames(plot_data), main=title,
                             ylab="Count",
                             ylim=c(0, max(colSums(plot_data, na.rm=TRUE))),
                             Legendtitle="Accumulation")
    }

    # If more samples, split samples up
    else{
        pages <- ceiling(ncol(plot_data) / 15)

        for(p in 1:pages){
            columns <- (((p - 1) * 15) + 1):(p * 15)
            columns <- columns[which(columns <= ncol(plot_data))]
            p <- Barplotsstacked(plot_data[, columns], AvgLine=FALSE,
                                 col=c("chocolate", "grey"),
                                 shownumbers_total=FALSE, shownumbers=TRUE,
                                 Legends=rownames(plot_data), main=title,
                                 ylab="Count",
                                 ylim=c(0, max(colSums(plot_data, na.rm=TRUE))),
                                 Legendtitle="Accumulation")
        }
    }

    # Plot per sample how many decoy feature quantifications were significant and
    # not significant
    loc <- which(grepl("_d$", rownames(temp_pval_quant)))
    plot_data <- rbind(colSums(temp_pval_quant[loc, ] < Quant_pVal_cut,
                               na.rm=TRUE),
                       colSums(temp_pval_quant[loc, ] > Quant_pVal_cut,
                               na.rm=TRUE))
    rownames(plot_data) <- c("certain", "uncertain")
    title <- "Certainty of ion accumulation for decoy features"

    # Only 15 samples can be plotted in one plot
    if(ncol(plot_data) <= 15){
        p <- Barplotsstacked(plot_data, AvgLine=FALSE, col=c("chocolate",
                                                             "grey"),
                             shownumbers_total=FALSE, shownumbers=TRUE,
                             Legends=rownames(plot_data), main=title,
                             ylab="Count",
                             ylim=c(0, max(colSums(plot_data, na.rm=TRUE))),
                             Legendtitle="Accumulation")
    }

    # If more samples, split samples up
    else{
        pages <- ceiling(ncol(plot_data) / 15)

        for(p in 1:pages){
            columns <- (((p - 1) * 15) + 1):(p * 15)
            columns <- columns[which(columns <= ncol(plot_data))]
            p <- Barplotsstacked(plot_data[, columns], AvgLine=FALSE,
                                 col=c("chocolate", "grey"),
                                 shownumbers_total=FALSE, shownumbers=TRUE,
                                 Legends=rownames(plot_data), main=title,
                                 ylab="Count",
                                 ylim=c(0, max(colSums(plot_data, na.rm=TRUE))),
                                 Legendtitle="Accumulation")
        }
    }


    #Boxplot per sample ratio of feature signal/background
    features_background_intensity <- log10(10^feat_w_backgr_int
                                           - 10^features_intensity)
    S2B <- base::log2(10^features_intensity / 10^features_background_intensity)
    S2B[is.infinite(as.matrix(S2B))] <- NA
    S2B <- base::as.data.frame(S2B)
    rownames(S2B) <- features$Feature_name

    lst <- list(Target_Decoy=ifelse(grepl("_d", rownames(S2B)), "decoy",
                                    "target"), S2B=S2B)
    QC_data[["Signal_to_background_target_decoy"]] <- lst
    loc <- which(!grepl("_d",rownames(S2B)))
    title <- "Signal to background ratio per feature quantification"

    # Only 15 samples can be plotted in one plot
    if(ncol(S2B) <= 15){
        graphics::par(mar=c(10, 4, 4, 4))
        graphics::boxplot(S2B[loc, ], outline=FALSE, main=title,
                          ylab="Signal/Background, log2", names=colnames(S2B),
                          las=2)
        graphics::abline(h=0, lty=2, col="red")

        for(co in 1:ncol(S2B)){
            med <- stats::median(S2B[loc, co], na.rm=TRUE)
            stxt <- ((graphics::par("usr")[4] - graphics::par("usr")[3]) * 0.05)
            graphics::text(co, med + stxt, round(med, digits=1))
        }
    }

    else{
        graphics::par(mar=c(10, 4, 4, 4))
        pages <- ceiling(ncol(S2B) / 15)

        for(p in 1:pages){
            columns <- (((p - 1) * 15) + 1):(p * 15)
            columns <- columns[which(columns <= ncol(S2B))]
            graphics::boxplot(S2B[loc, columns], outline=FALSE, main=title,
                              ylab="Signal/Background, log2",
                              names=colnames(S2B)[columns], las=2)
            graphics::abline(h=0, lty=2, col="red")

            for(co in 1:length(columns)){
                med <- stats::median(S2B[loc, columns[co]], na.rm=TRUE)
                stxt <- ((graphics::par("usr")[4] - graphics::par("usr")[3])
                        　* 0.05)
                graphics::text(co, med + stxt, round(med, digits=1))
            }
        }
    }

    pval_sig_w_bckgrnd_quant <- base::as.data.frame(pval_sig_w_bckgrnd_quant)
    features_intensity <- as.matrix(features_intensity)
    feat_w_backgr_int <- as.matrix(feat_w_backgr_int)
    features_intensity[is.infinite(features_intensity)] <- NA
    feat_w_backgr_int[is.infinite(feat_w_backgr_int)] <- NA
    features_intensity <- base::as.data.frame(features_intensity)
    feat_w_backgr_int <- base::as.data.frame(feat_w_backgr_int)

    # Calculate some scores which give insight in confidence of quantifications
    # (selection of right peaks)
    # Arrange selected peaks results in a long table
    peaks <- base::data.frame(
        RT=unlist(sapply(peak_quant, function(x){
            x$Peak_rt_with_background))
        },
        mz=unlist(sapply(peak_quant, function(x){
            x$Peak_mz_with_background))
        },
        known=unlist(sapply(peak_quant, function(x){
            x$correct_peak_w_background))
        },
        ion_count=unlist(sapply(peak_quant, function(x){
            x$Icount_feat_w_bkgr_int))
        },
        intensity=unlist(sapply(peak_quant, function(x){
            x$feat_w_backgr_int)
        })
    )
    lens <- length(samples)
    peaks$sample <- rep(sort(rep(samples, nrow(features))), length(peak_quant))
    peaks$feature <- rep(features$Feature_name, length(peak_quant) * lens)
    peaks$peak <- sort(rep(1:length(peak_quant), nrow(features) * lens))
    peaks$known[is.na(peaks$known)] <- 0

    # Add information which peak was finally used
    peaks$selected <- NA
    for(s in samples){
        selection_sample <- which(peaks$sample == s)
        matching_indices <- match(peaks$feature[selection_sample],
                                  rownames(peak_selected))
        peaks$selected[selection_sample] <- peak_selected[matching_indices, s]
    }

    # Add corrected RT and mz per entry (eliminate expected deviations between
    # samples)
    peaks$RT_correct <- NA
    peaks$mz_correct <- NA
    peaks$RT_correct <- as.numeric(peaks$RT_correct)
    peaks$mz_correct <- as.numeric(peaks$mz_correct)

    for(s in 1:length(samples)){
        aux <- base::paste("_calibration\\.", base::gsub("-", ".", samples[s]),
                           "$", sep="")
        corrections <- features[, c(1, which(grepl(aux, colnames(features))))]
        selection_peaks <- which(peaks$sample == samples[s])
        loc <- match(peaks$feature[selection_peaks], corrections$Feature_name)
        RT_corrected <- peaks$RT[selection_peaks] - corrections[loc, 2]
        mz_corrected <- peaks$mz[selection_peaks] - corrections[loc, 3]
        data.table::set(peaks, as.integer(selection_peaks), c(10L, 11L),
                        value=list(RT_corrected, mz_corrected))
    }


    # Calculate general peak variability score using delta_rt and delta_mz as sd
    # for zscoring. Feature with significant score here should be potentially
    # excluded
    Variability_alignment_scoring <- function(features, peaks, delta_rt,
                                              delta_mz,
                                              align_var_score_thr=0.05,
                                              plot=TRUE){
        # TODO: Remove (already declared in top level)
        # Convert z scores (vector or matrix) to pvalues
        zscore_to_pval <- function(z, alternative="greater"){
            if(alternative == "greater"){
                pval <- stats::pnorm(-z)
            }

            else if(alternative == "two.sided"){
                pval <- stats::pnorm(-abs(z))
                pval <- 2 * pval
            }

            else if(alternative == "less"){
                pval <- stats::pnorm(z)
            }

            return(pval)
        }

        # Select peaks which were selected by peak selection
        temp <- peaks[which(peaks$peak == peaks$selected), ]

        # Calculate interquartile range per feature over samples for RT and mz
        temp_RT <- stats::aggregate(temp$RT, by=list(Feature=temp$feature),
                                    FUN=matrixStats::iqr, na.rm=TRUE)
        temp_mz <- stats::aggregate(temp$mz, by=list(Feature=temp$feature),
                                    FUN=matrixStats::iqr, na.rm=TRUE)

        # Z score using delta RT and delta mz as sd
        temp_RT$x <- (temp_RT$x - mean(temp_RT$x, na.rm=TRUE)) / delta_rt
        temp_mz$x <- (temp_mz$x - mean(temp_mz$x, na.rm=TRUE)) / delta_mz

        # Features with a negative z-score represent features with much lower
        # deviation between samples than expected
        # Features with a positive z-score represent features with much higher
        # deviation between samples than expected --> have to be inspected and
        # eventually removed
        align_var_score <- base::data.frame(
            RT_variability_pval=(zscore_to_pval(temp_RT$x)),
            mz_variability_pval=(zscore_to_pval(temp_mz$x))
        )
        loc <- match(features$Feature_name, temp_RT$Feature)
        align_var_score <- align_var_score[loc, ]
        rownames(align_var_score) <- features$Feature_name
        align_var_score <- base::as.data.frame(align_var_score)
        rmin <- matrixStats::rowMins(as.matrix(align_var_score), na.rm=TRUE)
        align_var_score <- cbind(align_var_score, rmin)
        colnames(align_var_score)[3] <- "combined_variability_pval"

        # Add some plots showing number of features with general high variability
        # between samples
        if(plot == TRUE){
            loc <- which(!grepl("_d|_i", rownames(align_var_score)))
            temp <- matrixStats::rowMins(as.matrix(align_var_score[loc, ]),
                                         na.rm=TRUE)
            x1 <- which(temp > align_var_score_thr)
            x2 <- which(temp < align_var_score_thr)
            Variability_scores <- cbind(length(x1), length(x2))
            colnames(Variability_scores) <- c("normal", "high")
            title <- "Variability in peak selection between samples"
            p <- BarplotsSBS(Variability_scores, AvgLine=FALSE,
                             col=c("chocolate", "grey"), shownumbers=TRUE,
                             main=title, ylab="Count")
        }

        return(align_var_score)
    }

    align_var_score <- Variability_alignment_scoring(features, peaks, delta_rt,
                                                     delta_mz,
                                                     align_var_score_thr)

    # Calculate standardized (z score) deviation of peak RT/mz from mean over
    # all samples with used sd defined by delta_RT and delta_mz. Z scores are
    # then translated into two sided pvalues and per feature/sample the
    # minimum of pval of RT or mz is reported
    Alignment_scoring <- function(peaks, features_select, samples, sd_RT=0.5,
                                  sd_mz=0.001, corrected_alignment=TRUE,
                                  align_score_thr=0.05, plot=TRUE){
        # Determine for every sample if detected peak RT or mz is significantly
        # deviating from all other determined peak RT or mz
        # Extract only finally used peaks
        temp_peaks <- peaks[which(peaks$peak == peaks$selected), ]
        # Order according to feature name
        temp_peaks <- temp_peaks[order(temp_peaks$feature, temp_peaks$sample), ]

        # Prepare matrices which should store RT and mz of selected peaks per
        # sample and feature
        if(corrected_alignment == TRUE){
            loc <- c("feature", "sample", "RT_correct")
        }

        if(corrected_alignment == FALSE){
            loc <- c("feature", "sample", "RT")
        }

        sel_peaks_RT <- stats::reshape(temp_peaks[, loc], idvar="feature",
                                            timevar="sample", direction="wide")

        rownames(sel_peaks_RT) <- sel_peaks_RT$feature
        sel_peaks_RT <- sel_peaks_RT[, -1]

        if(corrected_alignment == TRUE){
            loc <- c("feature", "sample", "mz_correct")
        }

        if(corrected_alignment == FALSE){
            loc <- c("feature", "sample", "mz")
        }

        sel_peaks_mz <- stats::reshape(temp_peaks[, loc], idvar="feature",
                                       timevar="sample", direction="wide")

        rownames(sel_peaks_mz) <- sel_peaks_mz$feature
        sel_peaks_mz <- sel_peaks_mz[, -1]

        # Calculate mean and sd of RT and mz per feature
        mn <- rowMeans(as.matrix(sel_peaks_RT), na.rm=TRUE)
        sd <- matrixStats::rowSds(as.matrix(sel_peaks_RT, na.rm=TRUE))
        sel_peaks_dists_RT <- base::data.frame(mean=, sd=sd)
        mn <- rowMeans(as.matrix(sel_peaks_mz), na.rm=TRUE)
        sd <- matrixStats::rowSds(as.matrix(sel_peaks_mz, na.rm=TRUE))
        sel_peaks_dists_mz <- base::data.frame(mean=mn, sd=sd)

        # Convert observed RTs and mzs per feature and sample into z-scores
        sel_peaks_RT_zscore <- (sel_peaks_RT - sel_peaks_dists_RT$mean) / sd_RT
        sel_peaks_mz_zscore <- (sel_peaks_mz - sel_peaks_dists_mz$mean) / sd_mz

        # Convert zscores into pvalues
        zscore_to_pval <- function(z){
            pval <- stats::pnorm(-abs(z))
            pval <- 2 * pval

            return(pval)
        }

        sel_peaks_RT_pval <- base::as.data.frame(apply(sel_peaks_RT_zscore, 1:2,
                                                       zscore_to_pval))
        sel_peaks_mz_pval <- base::as.data.frame(apply(sel_peaks_mz_zscore, 1:2,
                                                       zscore_to_pval))
        cname1 <- colnames(sel_peaks_RT_pval)
        cname2 <- colnames(sel_peaks_mz_pval)

        if(corrected_alignment == TRUE){
            colnames(sel_peaks_RT_pval) <- base::gsub("RT_correct.", "", cname1)
            colnames(sel_peaks_mz_pval) <- base::gsub("mz_correct.", "", cname2)
        }

        if(corrected_alignment == FALSE){
            cnames <-
            colnames(sel_peaks_RT_pval) <- base::gsub("RT.", "", cname1)
            colnames(sel_peaks_mz_pval) <- base::gsub("mz.", "", cname2)
        }

        # Reorder rows to match to ordering in features_select
        loc1 <- match(features_select$Feature_name, rownames(sel_peaks_RT_pval))
        loc2 <- match(features_select$Feature_name, rownames(sel_peaks_mz_pval))
        sel_peaks_RT_pval <- sel_peaks_RT_pval[loc1, ]
        sel_peaks_mz_pval <- sel_peaks_mz_pval[loc2, ]

        # Summarize alignment score per feature and sample to single minimal
        # pval (either RT or mz)
        alignment_score_peaks <- pmin(sel_peaks_RT_pval, sel_peaks_mz_pval)

        # Add some plots about for how many features we are certain/uncertain
        # about the peak selection
        if(plot == TRUE){
            loc <- which(!grepl("_d|_d_i", rownames(alignment_score_peaks)))
            temp_score <- alignment_score_peaks[loc, ]
            temp_score <- rbind(colSums(temp_score > align_score_thr,
                                        na.rm=TRUE),
                                colSums(temp_score < align_score_thr,
                                        na.rm=TRUE))
            rownames(temp_score) <- c("certain", "uncertain")

            main <- ifelse(corrected_alignment == TRUE,
                           "Certainty of peak selection (corrected)",
                           "Certainty of peak selection (raw)")

            # Only 15 samples can be plotted in one plot
            if(ncol(temp_score) <= 15){
                p <- Barplotsstacked(temp_score, AvgLine=FALSE,
                                     col=c("chocolate", "grey"),
                                     shownumbers_total=FALSE, shownumbers=TRUE,
                                     Legends=rownames(temp_score), main=main,
                                     ylab="Count", Legendtitle="Peak selection",
                                     ylim=c(0, max(colSums(temp_score,
                                                           na.rm=TRUE))))
            }

            # If more samples, split samples up
            else{
                pages <- ceiling(ncol(temp_score) / 15)
                for(p in 1:pages){
                    columns <- (((p - 1) * 15) + 1):(p * 15)
                    columns <- columns[which(columns <= ncol(temp_score))]
                    p <- Barplotsstacked(temp_score[, columns], AvgLine=FALSE,
                                         col=c("chocolate", "grey"),
                                         shownumbers_total=FALSE,
                                         shownumbers=TRUE,
                                         Legends=rownames(temp_score),
                                         main=main, ylab="Count",
                                         Legendtitle="Peak selection",
                                         ylim=c(0, max(colSums(temp_score,
                                                               na.rm=TRUE))))
                }
            }
        }

        return(alignment_score_peaks)
    }

    # Perform alignment scoring based on corrected mz and RT per feature/sample
    align_scores_peaks_correct <- Alignment_scoring(
        peaks,
        features,
        samples,
        sd_RT=delta_rt,
        sd_mz=delta_mz,
        corrected_alignment=TRUE,
        align_score_thr=align_score_thr
    )
    # Perform alignment scoring based on raw mz and RT per feature/sample to
    # detect unexpected complete outliers
    align_scores_peaks_raw <- Alignment_scoring(peaks, features, samples,
                                                sd_RT=delta_rt, sd_mz=delta_mz,
                                                corrected_alignment=FALSE,
                                                align_score_thr=align_score_thr)

    # Calculate score based on cooccurence of isotope peaks
    c1 <- grepl("_i", features$Feature_name)
    c2 <- !grepl("_d_i", features$Feature_name)

    if(any(c1 & c2)){
        # XXX: Defining a function inside a conditional?
        Mono_iso_alignment_scoring <- function(features, samples, peaks,
                                               pval_sig_w_bckgrnd_quant,
                                               delta_rt, delta_mz,
                                               mono_iso_alignment_cutoff,
                                               plot=TRUE){
            #prepare dataframe into which all results are saved
            mat <- matrix(nrow=nrow(features), ncol=length(samples))
            mono_iso_alignment_summary <- base::as.data.frame(mat)
            rownames(mono_iso_alignment_summary) <- features$Feature_name

            for(c in 1:ncol(mono_iso_alignment_summary)){
                mias <- mono_iso_alignment_summary[, c]
                mono_iso_alignment_summary[, c] <- as.numeric(mias)
            }

            if(length(which(peaks$known == 1)) >= 10){
                # XXX: Again?
                zscore_to_pval <- function(z, alternative="greater"){
                    if(alternative == "greater"){
                        pval <- stats::pnorm(-z)
                    }

                    else if(alternative == "two.sided"){
                        pval <- stats::pnorm(-abs(z))
                        pval <- 2 * pval
                    }

                    else if(alternative=="less"){
                        pval <- stats::pnorm(z)
                    }

                    return(pval)
                }

                # Estimate expected distributions of deviations in RT and mz
                # between M and M+1 peaks
                mean_RT_distribution <- vector("numeric", length(samples))
                sd_RT_distribution <- vector("numeric", length(samples))
                mean_mz_distribution <- vector("numeric", length(samples))
                sd_mz_distribution <- vector("numeric", length(samples))

                for(s in 1:length(samples)){
                    temp <- peaks[which(peaks$sample == samples[s]), ]
                    # Add quantification scores
                    loc <- match(temp$feature,
                                 rownames(pval_sig_w_bckgrnd_quant))
                    temp$quant_score <- pval_sig_w_bckgrnd_quant[loc, s]
                    # Add charge state
                    loc <- match(features$Feature_name, temp$feature)
                    temp$charge <- features$Charge[loc]
                    # Selected peaks RT and mz
                    temp_selected <- temp[which(temp$peak == temp$selected), ]
                    # Split features into mono and +1 isotopes
                    feat <- temp_selected$feature
                    temp_selected_mono <- temp_selected[which(!grepl("_i",
                                                                     feat)), ]
                    temp_selected_iso <- temp_selected[which(grepl("_i",
                                                                   feat)), ]
                    feat <- temp_selected_iso$feature
                    temp_selected_iso$feature <- base::gsub("_i", "", feat)
                    # Combine both and only keep features for which we also have
                    # +1 isotope features
                    combined <- dplyr::inner_join(temp_selected_mono,
                                                  temp_selected_iso,
                                                  by="feature")
                    # Determine expected variation between mono and +1 isotopes
                    # based on true identifications and significantly quantified
                    # mono and +1 isotope features
                    true <- combined[which(combined$known.x == 1
                                           & combined$quant_score.x < 0.05
                                           & combined$quant_score.y < 0.05), ]
                    true$delta_rt_iso_mono <- true$RT.y - true$RT.x
                    aux <- ((true$mz.y * true$charge.y) - 1.002054)
                    true$delta_mz_iso_mono <- (aux / true$charge.y) - true$mz.x

                    mean_RT_distribution[s] <- mean(true$delta_rt_iso_mono,
                                                    na.rm=TRUE)
                    sd_RT_distribution[s] <- stats::sd(true$delta_rt_iso_mono,
                                                       na.rm=TRUE)
                    mean_mz_distribution[s] <- mean(true$delta_mz_iso_mono,
                                                    na.rm=TRUE)
                    sd_mz_distribution[s] <- stats::sd(true$delta_mz_iso_mono,
                                                       na.rm=TRUE)
                }

                # Estimate normal distribution of deviations in RT and mz
                # between mono and +1 isotope peak
                x <- seq(-delta_rt * 3, delta_rt * 3, length=1000)
                y <- stats::dnorm(x, mean=mean(mean_RT_distribution,
                                               na.rm=TRUE),
                                  sd=mean(sd_RT_distribution, na.rm=TRUE))
                RT_deviation_norm <- base::data.frame(x=x, y=y)
                x <- seq(-delta_mz * 3, delta_mz * 3, length=1000)
                y <- stats::dnorm(x, mean=mean(mean_mz_distribution,
                                               na.rm=TRUE),
                                  sd=mean(sd_mz_distribution, na.rm=TRUE))
                mz_deviation_norm <- base::data.frame(x=x, y=y)

                title <- "RT-deviation M vs M+1 peaks - true quantifications"
                graphics::plot(RT_deviation_norm$x, RT_deviation_norm$y,
                               type="l", xlab="RT(M+1) - RT(M)", main=title,
                               ylab="Density")
                graphics::abline(v=mean(true$delta_rt_iso_mono), lty=2)
                title <- "mz-deviation M vs M+1 peaks - true quantifications"
                graphics::plot(mz_deviation_norm$x, mz_deviation_norm$y,
                               type="l", xlab="mz(M+1) - mz(M)", main=title,
                               ylab="Density")
                graphics::abline(v=mean(true$delta_mz_iso_mono), lty=2)

                # Now calculate for every available feature for which an isotope
                # feature is available how well both peaks are aligned
                mean_RT_distribution_mean <- mean(mean_RT_distribution,
                                                  na.rm=TRUE)
                sd_RT_distribution_mean <- mean(sd_RT_distribution, na.rm=TRUE)
                mean_mz_distribution_mean <- mean(mean_mz_distribution,
                                                  na.rm=TRUE)
                sd_mz_distribution_mean <- mean(sd_mz_distribution, na.rm=TRUE)

                for(s in 1:length(samples)){
                    temp <- peaks[which(peaks$sample == samples[s]), ]
                    # Add quantification scores
                    loc <- match(temp$feature,
                                 rownames(pval_sig_w_bckgrnd_quant))
                    temp$quant_score <- pval_sig_w_bckgrnd_quant[loc, s]
                    # Add charge state
                    temp$charge <- features$Charge[match(features$Feature_name,
                                                         temp$feature)]
                    # Selected peaks RT and mz
                    temp_selected <- temp[which(temp$peak == temp$selected), ]
                    # Split features into mono and +1 isotopes
                    feat <- temp_selected$feature
                    temp_selected_mono <- temp_selected[which(!grepl("_i",
                                                                     feat)), ]
                    temp_selected_iso <- temp_selected[which(grepl("_i",
                                                                   feat)), ]
                    feat <- temp_selected_iso$feature
                    temp_selected_iso$feature <- base::gsub("_i", "", feat)
                    # Combine both and only keep features for which we also have
                    # +1 isotope features
                    combined <- dplyr::inner_join(temp_selected_mono,
                                                  temp_selected_iso,
                                                  by="feature")
                    combined$delta_rt_iso_mono <- combined$RT.y - combined$RT.x
                    aux <- (((combined$mz.y * combined$charge.y) - 1.002054)
                            / combined$charge.y)
                    combined$delta_mz_iso_mono <- aux - combined$mz.x
                    # Convert to z-score
                    aux <- ((combined$delta_rt_iso_mono
                             - mean_RT_distribution_mean)
                            /sd_RT_distribution_mean)
                    combined$z_delta_rt_iso_mono <- aux
                    aux <- ((combined$delta_mz_iso_mono
                             -mean_mz_distribution_mean)
                            /sd_mz_distribution_mean)
                    combined$z_delta_mz_iso_mono <- aux
                    # Convert to pvalues
                    rt <- combined$z_delta_rt_iso_mono
                    mz <- combined$z_delta_mz_iso_mono
                    combined$RT_deviation_pval <- zscore_to_pval(rt,
                                                                 "two.sided")
                    combined$mz_deviation_pval <- zscore_to_pval(mz,
                                                                 "two.sided")
                    # Store results
                    mat <- as.matrix(combined[, c("RT_deviation_pval",
                                                  "mz_deviation_pval")])
                    combined$score <- matrixStats::rowMins(mat, na.rm=TRUE)
                    aux <- base::gsub("_i", "",
                                      rownames(mono_iso_alignment_summary))
                    loc <- match(aux, combined$feature)
                    mono_iso_alignment_summary[, s] <- combined$score[loc]
                    colnames(mono_iso_alignment_summary)[s] <- samples[s]
                }

                loc <- is.infinite(as.matrix(mono_iso_alignment_summary))
                mono_iso_alignment_summary[loc] <- NA

                # Add some plots about for how many features we are certain
                # about peak selection based on detected deviations to isotope
                # peak selection
                if(plot == TRUE){
                    temp_alignmentscores <- mono_iso_alignment_summary
                    loc <- which(!grepl("_d|_i|_d_i",
                                        rownames(temp_alignmentscores)))
                    temp_alignmentscores <- temp_alignmentscores[loc, ]
                    aux <- temp_alignmentscores > mono_iso_alignment_cutoff
                    cs1 <- colSums(aux, na.rm=TRUE)
                    aux <- temp_alignmentscores < mono_iso_alignment_cutoff
                    cs2 <- colSums(aux, na.rm=TRUE)
                    RT_alignment_scores <- rbind(cs1, cs2)
                    rownames(RT_alignment_scores) <- c("certain", "uncertain")

                    # Only 15 samples can be plotted in one plot
                    if(ncol(RT_alignment_scores) <= 15){
                        title <- paste0("Certainty of peak selection based on ",
                                        "mono/+1 isotope peak selection")
                        mx <- max(colSums(RT_alignment_scores, na.rm=TRUE))
                        rnames <- rownames(RT_alignment_scores)
                        p <- Barplotsstacked(RT_alignment_scores, AvgLine=FALSE,
                                             col=c("chocolate", "grey"),
                                             shownumbers_total=FALSE,
                                             shownumbers=TRUE,
                                             Legends=rnames,
                                             main=title, ylab="Count",
                                             ylim=c(0, mx),
                                             Legendtitle="Peak selection")
                    }

                    # If more samples, split samples up
                    else{
                        pages <- ceiling(ncol(RT_alignment_scores) / 15)

                        for(p in 1:pages){
                            columns <- (((p - 1) * 15) + 1):(p * 15)
                            loc <- which(columns <= ncol(RT_alignment_scores))
                            columns <- columns[loc]
                            legend <- rownames(RT_alignment_scores)
                            title <- paste0("Certainty of peak selection ",
                                            "based on mono/+1 isotope peak ",
                                            "selection")
                            mx <- max(colSums(RT_alignment_scores, na.rm=TRUE))
                            p <- Barplotsstacked(RT_alignment_scores[, columns],
                                                 AvgLine=FALSE,
                                                 col=c("chocolate", "grey"),
                                                 shownumbers_total=FALSE,
                                                 shownumbers=TRUE,
                                                 Legends=legend, main=title,
                                                 ylab="Count", ylim=c(0, mx),
                                                 Legendtitle="Peak selection")
                        }
                    }
                }
            }

            else{
                print(paste0("Warning: The true peak for too few features is ",
                             "known. Skipping +1-isotope alignment scoring."))
            }

            return(mono_iso_alignment_summary)
        }

        mono_iso_alignment_summary <- Mono_iso_alignment_scoring(
            features,
            samples,
            peaks,
            pval_sig_w_bckgrnd_quant,
            delta_rt,
            delta_mz,
            mono_iso_alignment_cutoff
        )
    }

    else{
        mono_iso_alignment_summary <- NA
    }

    FDR_peak_selection <- Peak_selection_FDR(num_features=500, features,
                                             samples, peaks, path_to_features,
                                             peak_quant, feat_w_backgr_int,
                                             peak_selected, delta_mz, delta_rt,
                                             peak_min_ion_count,
                                             num_peaks_store, align_score_thr,
                                             n_cores, peak_decision,
                                             Alignment_scoring, plot=TRUE)

    QC_data$FDR_peak_selection <- FDR_peak_selection

    if(any(!is.na(FDR_peak_selection$Large_Intensity_delta_FDR))){
        for(i in 1:length(FDR_peak_selection$Total_FDR)){
            rnd <- round(FDR_peak_selection$Large_Intensity_delta_FDR[i],
                         digits=1)
            print(base::paste(samples[i], " - FDR of peak selection: ", rnd,
                              " %", sep=""))
        }
    }

    grDevices::dev.off()

    print(paste0(Sys.time(),
                 " Finished scoring of alignment and quantification"))

    # Evaluate alignment performance
    # After peak alignment
    peaks_selected <- peaks[which(peaks$peak == peaks$selected), ]
    # Raw peaks
    peaks_raw <- peaks[which(peaks$peak == 6), ]
    select_only_known <- which(base::paste(peaks_raw$sample, peaks_raw$feature)
                               %in% base::paste(peaks_selected$sample,
                                                peaks_selected$feature))
    peaks_raw <- peaks_raw[select_only_known, ]

    # Deviation in raw data
    usecols <- c("RT", "mz", "RT_correct", "mz_correct")
    temp_mean_raw <- stats::aggregate(peaks_raw[, usecols],
                                      by=list(peaks_raw$feature), FUN=mean,
                                      na.rm=TRUE)
    temp_sd_raw <- stats::aggregate(peaks_raw[, usecols],
                                    by=list(peaks_raw$feature), FUN=stats::sd,
                                    na.rm=TRUE)
    # Deviation after alignment
    temp_mean_aligned <- stats::aggregate(peaks_selected[, usecols],
                                          by=list(peaks_selected$feature),
                                          FUN=mean, na.rm=TRUE)
    temp_sd_aligned <- stats::aggregate(peaks_selected[, usecols],
                                        by=list(peaks_selected$feature),
                                        FUN=stats::sd, na.rm=TRUE)
    # Exclude features where sd in uncorrected RT > 1 min
    sel <- which(temp_sd_raw$RT > 1)
    temp_mean_raw <- temp_mean_raw[-sel, ]
    temp_sd_raw <- temp_sd_raw[-sel, ]
    temp_mean_aligned <- temp_mean_aligned[-sel, ]
    temp_sd_aligned <- temp_sd_aligned[-sel, ]
    grDevices::pdf("Performance of feature alignment.pdf", useDingbats=FALSE)

    tryCatch(
        {
            graphics::boxplot(temp_sd_raw$RT, temp_sd_aligned$RT_correct,
                              outline=FALSE, main="Variability of feature RT",
                              names=c("Raw", "Aligned"), ylab="SD of RT [min]")
            graphics::boxplot(temp_sd_raw$mz, temp_sd_aligned$mz_correct,
                              outline=FALSE, main="Variability of feature mz",
                              names=c("Raw", "Aligned"), ylab="SD of mz [Da]")
        },
        error=function(cond){
            warning(paste("Error while plotting performance of feature ",
                          "alignment. Original message \n", cond))
        }
    )

    grDevices::dev.off()

    # Save temporary results
    print(paste0(Sys.time(),
                 " Save quantification results after alignment scoring"))
    temp_results <- list(features=features,
                         features_intensity=features_intensity,
                         Icount_feat_sample_mat=Icount_feat_sample_mat,
                         feat_w_backgr_int=feat_w_backgr_int,
                         Icount_feat_w_bkgr_int=Icount_feat_w_bkgr_int,
                         peak_selected=peak_selected, S2B=S2B,
                         pval_sig_w_bckgrnd_quant=pval_sig_w_bckgrnd_quant,
                         align_var_score=align_var_score,
                         align_scores_peaks_correct=align_scores_peaks_correct,
                         align_scores_peaks_raw=align_scores_peaks_raw,
                         mono_iso_alignment_summary=mono_iso_alignment_summary,
                         peaks=peaks, QC_data)
    setwd(base::paste(path_to_features, "/Temporary_files", sep=""))
    save(temp_results, file="Quantification_raw_results_with_scores.RData")

    crap <- gc(FALSE)

    sel <- which(!grepl("_d", features$Feature_name))

    if(length(selection) > 0){
        features <- features[sel, ]
        Icount_feat_sample_mat <- Icount_feat_sample_mat[sel, ]
        feat_w_backgr_int <- feat_w_backgr_int[sel, ]
        Icount_feat_w_bkgr_int <- Icount_feat_w_bkgr_int[sel, ]
        peak_selected <- peak_selected[sel, ]
        S2B <- S2B[sel, ]
        pval_sig_w_bckgrnd_quant <- pval_sig_w_bckgrnd_quant[sel, ]
        align_var_score <- align_var_score[sel, ]
        align_scores_peaks_correct <- align_scores_peaks_correct[sel, ]
        align_scores_peaks_raw <- align_scores_peaks_raw[sel, ]

        if(any(!is.na(mono_iso_alignment_summary))){
            mono_iso_alignment_summary <- mono_iso_alignment_summary[sel, ]
        }
    }

    # Remove isotope features which never show significant quantification
    sel <- which(grepl("_i", features$Feature_name)
                 & matrixStats::rowMins(as.matrix(pval_sig_w_bckgrnd_quant),
                                        na.rm=TRUE) > Quant_pVal_cut)

    total_length <- length(which(grepl("_i", features$Feature_name)))

    if(length(selection) > 0){
        features <- features[-sel, ]
        Icount_feat_sample_mat <- Icount_feat_sample_mat[-sel, ]
        feat_w_backgr_int <- feat_w_backgr_int[-sel, ]
        Icount_feat_w_bkgr_int <- Icount_feat_w_bkgr_int[-sel, ]
        peak_selected <- peak_selected[-sel, ]
        S2B <- S2B[-sel, ]
        pval_sig_w_bckgrnd_quant <- pval_sig_w_bckgrnd_quant[-sel, ]
        align_var_score <- align_var_score[-sel, ]
        align_scores_peaks_correct <- align_scores_peaks_correct[-sel, ]
        align_scores_peaks_raw <- align_scores_peaks_raw[-sel, ]

        if(any(!is.na(mono_iso_alignment_summary))){
            mono_iso_alignment_summary <- mono_iso_alignment_summary[-sel, ]
        }

        num <- round(length(selection) / total_length * 100, digits=1)
        print(base::paste(Sys.time(), "Removed ", length(selection),
                          " isotope features (", num, " %) as they dont show ",
                          "significant ion accumulation in any sample.",
                          sep=""))
    }

    # Impute missing values based on generalized additive model
    impute_feature_level_quant <- function(data, features,
                                           backgr_inty_GAM_tab_x_feat_sum){
        quant_data <- data.table::copy(data)
        backgr_inty_GAM_tab_x_feat <- backgr_inty_GAM_tab_x_feat_sum

        for(c in 1:ncol(quant_data)){
            missing_vals_rows <- which(is.na(quant_data[, c]))
            loc <- match(rownames(quant_data)[missing_vals_rows],
                         rownames(backgr_inty_GAM_tab_x_feat))
            impute_vals <- backgr_inty_GAM_tab_x_feat[loc, c]
            impute_vals <- log10(2^impute_vals)
            data.table::set(quant_data, i=as.integer(missing_vals_rows),
                            j=as.integer(c), value=impute_vals)
        }

        return(quant_data)
    }

    gtab <- QC_data$Decoy_feature_parameters$backgr_inty_GAM_tab_x_feat_sum
    feat_w_backgr_int_imputed <- impute_feature_level_quant(
        data=feat_w_backgr_int,
        features=features,
        backgr_inty_GAM_tab_x_feat_sum=gtab
    )

    # Now perform filtering of quantifications based on scores
    # remove all features with too high variability
    loc <- which(matrixStats::rowMins(as.matrix(align_var_score), na.rm=TRUE)
                 < align_var_score_thr)
    selection <- rownames(align_var_score)[loc]

    if(length(selection) > 0){
        print(base::paste(Sys.time(), "Removed ",
                          length(which(!grepl("_i", selection))), " (",
                          length(which(grepl("_i",selection))), " isotope) ",
                          "features from quantification results due to too ",
                          "high variability in alignment between samples.",
                          sep=""))
        sel <- match(selection, features$Feature_name)

        features <- features[-sel, ]
        Icount_feat_sample_mat <- Icount_feat_sample_mat[-sel, ]
        feat_w_backgr_int <- feat_w_backgr_int[-sel, ]
        feat_w_backgr_int_imputed <- feat_w_backgr_int_imputed[-sel, ]
        Icount_feat_w_bkgr_int <- Icount_feat_w_bkgr_int[-sel, ]
        peak_selected <- peak_selected[-sel, ]
        S2B <- S2B[-sel, ]
        pval_sig_w_bckgrnd_quant <- pval_sig_w_bckgrnd_quant[-sel, ]
        align_var_score <- align_var_score[-sel, ]
        align_scores_peaks_correct <- align_scores_peaks_correct[-sel, ]
        align_scores_peaks_raw <- align_scores_peaks_raw[-sel, ]

        if(any(!is.na(mono_iso_alignment_summary))){
            mono_iso_alignment_summary <- mono_iso_alignment_summary[-sel, ]
        }
    }

    # Remove all quantifications per sample with uncertain peak selection
    for(s in 1:length(samples)){
        loc <- which(grepl(samples[s], colnames(align_scores_peaks_correct)))
        temp_scores <- align_scores_peaks_correct[, loc, drop=FALSE]
        temp_scores_raw <- align_scores_peaks_raw[, loc, drop=FALSE]
        sel <- which(temp_scores < align_score_thr)
        names(sel) <- rownames(temp_scores)[sel]
        feat_w_backgr_int[sel, s] <- NA
        feat_w_backgr_int_imputed[sel, s] <- NA

        num <- round(length(sel) / nrow(align_scores_peaks_correct) * 100,
                     digits=1)
        print(base::paste(samples[s], ": Removed ",
                          length(which(!grepl("_i", names(sel)))), " (",
                          length(which(grepl("_i", names(sel)))),
                          " isotope) quantifications (", num,
                          " %) due to uncertain peak selection", sep=""))
    }

    # Remove isotope quantifications for uncertain mono-+1-isos but keep mono
    # quantification if otherwise certain of peak selection
    if(any(!is.na(mono_iso_alignment_summary))){
        for(s in 1:length(samples)){
            if(length(which(peaks$known == 1)) >= 10){
                loc <- which(grepl(samples[s],
                                   colnames(mono_iso_alignment_summary)))
                temp_scores <- mono_iso_alignment_summary[, loc]
                sel <- which(temp_scores < mono_iso_alignment_cutoff
                             & grepl("_i",
                                     rownames(mono_iso_alignment_summary)))
                names(selection) <- rownames(mono_iso_alignment_summary)[sel]

                # How many of these still have a quantification (not already
                # removed by previous filtering step)
                temp_quant <- feat_w_backgr_int[sel, s]
                feat_w_backgr_int[sel, s] <- NA
                feat_w_backgr_int_imputed[sel, s] <- NA
                num <- round(length(which(!is.na(temp_quant)))
                             / nrow(align_scores_peaks_correct) * 100, digits=1)
                print(base::paste(samples[s], ": Removed ",
                                  length(which(!is.na(temp_quant))),
                                  " isotope quantifications (", num,
                                  " %) due to discrepancy in peak selection ",
                                  "between mono-/+1-isotope features", sep=""))
            }

            # Not enough known peaks ... remove all isotopes
            else{
                selection <- which(grepl("_i", rownames(feat_w_backgr_int)))
                feat_w_backgr_int[sel, s] <- NA
                feat_w_backgr_int_imputed[sel, s] <- NA

                print(base::paste(samples[s], ": Removed all isotope ",
                                  "quantifications due to too few known true ",
                                  "peaks.", sep=""))
            }
        }
    }

    # Calculate fraction of missing values before and after imputation
    total <- nrow(feat_w_backgr_int) * ncol(feat_w_backgr_int)
    missing <- length(which(is.na(as.numeric(as.matrix(feat_w_backgr_int)))))
    print(base::paste(Sys.time(), "Quantification without imputation: ",
                      round(missing / total * 100, digits=1),
                      " % missing values", sep=""))
    aux <- as.numeric(as.matrix(feat_w_backgr_int_imputed))
    missing <- length(which(is.na(aux)))
    print(base::paste(Sys.time(), "Quantification with imputation: ",
                      round(missing / total * 100, digits=1),
                      " % missing values", sep=""))

    # Store results after filtering
    print(paste0(Sys.time(), " Save quantification results after alignment ",
                 "scoring and filter"))

    temp_results <- list(features=features,
                         features_intensity=features_intensity,
                         Icount_feat_sample_mat=Icount_feat_sample_mat,
                         feat_w_backgr_int=feat_w_backgr_int,
                         feat_w_backgr_int_imputed=feat_w_backgr_int_imputed,
                         Icount_feat_w_bkgr_int=Icount_feat_w_bkgr_int,
                         peak_selected=peak_selected, S2B=S2B,
                         pval_sig_w_bckgrnd_quant=pval_sig_w_bckgrnd_quant,
                         align_var_score=align_var_score,
                         align_scores_peaks_correct=align_scores_peaks_correct,
                         align_scores_peaks_raw=align_scores_peaks_raw,
                         mono_iso_alignment_summary=mono_iso_alignment_summary,
                         QC_data)
    setwd(base::paste(path_to_features, "/Temporary_files", sep=""))
    save(temp_results,
         file="Quantification_raw_results_with_scores_filtered.RData")

    # Perform protein level quantification

    correct_intensities <- function(features, feat_sample_mat_requant,
                                    pval_quant, MaxQ_peptides_quant, main="",
                                    corr_factor=NA){
        # No correction factor specified so correction will be determined and
        # then applied
        if(is.na(corr_factor)){
            # Select monoisotopic features with significant quantification for
            # determination of general trends in difference between MaxQ and
            # Requant quantification
            select <- which(!grepl("_i|_d", features$Feature_name))
            feature_quant <- feat_sample_mat_requant[select, ]
            feature_quant_pval <- pval_quant[select, ]
            loc <- is.na(feature_quant_pval) | feature_quant_pval > 0.1
            feature_quant[loc] <- NA
            features_select <- features[select, ]
            # Determine deviations between Requant and MaxQ in abundance
            # estimations
            lst <- list(Sequence=features_select$Sequence)
            Requant_peptides_quant_seq <- stats::aggregate(10^feature_quant,
                                                           by=lst, FUN=sum)
            logs <- base::log2(Requant_peptides_quant_seq[, -1])
            Requant_peptides_quant_seq[, -1] <- logs

            # Bring MaxQ results table into same order as Requant output
            if(ncol(MaxQ_peptides_quant) != ncol(Requant_peptides_quant_seq)){
                cnam <- colnames(MaxQ_peptides_quant)
                colnames(MaxQ_peptides_quant) <- base::gsub("Intensity.", "",
                                                            cnam)
                ordering <- vector("numeric", ncol(MaxQ_peptides_quant))

                for(i in 1:ncol(MaxQ_peptides_quant)){
                    cnames <- colnames(MaxQ_peptides_quant)[i]
                    overlap <- grepl(base::paste(cnames, "$", sep=""),
                                     colnames(Requant_peptides_quant_seq))
                    ordering[i] <- ifelse(any(overlap), which(overlap == T), NA)
                }

                loc <- !is.na(ordering)
                MaxQ_peptides_quant <- MaxQ_peptides_quant[, which(loc)]
                ordering <- ordering[loc]
                MaxQ_peptides_quant <- MaxQ_peptides_quant[, ordering]
            }

            comb <- dplyr::full_join(MaxQ_peptides_quant,
                                     Requant_peptides_quant_seq, by="Sequence")

            ncolumns <- ncol(feat_sample_mat_requant)

            dat <- base::as.data.frame(matrix(ncol=2,
                                              nrow=nrow(comb) * ncolumns))
            dat[, 1] <- as.numeric(dat[, 1])
            dat[, 2] <- as.numeric(dat[, 2])

            for(i in 1:ncolumns){
                start_index <- (i - 1) * nrow(comb) + 1
                end_index <- i * nrow(comb)
                temp <- base::as.data.frame(cbind(comb[, i + 1],
                                                  comb[, i + 1 + ncolumns]))
                data.table::set(dat, as.integer(start_index:end_index),
                                j=as.integer(1:2), value=temp)
            }

            colnames(dat) <- c("MaxQ", "Requant")
            dat <- stats::na.omit(dat)

            if(main != ""){
                main <- base::paste(" ", main, sep="")
            }

            grDevices::pdf(base::paste("Correct feature abundance estimations",
                                       main, ".pdf", sep=""))
            title <- "Correlation of MaxQ vs Requant on peptide level"
            graphics::smoothScatter(dat[, 1], dat[, 2], ylab="Requant, log2",
                                    xlab="MaxQ, log2", main=title)
            graphics::abline(a=0, b=1)
            nr <- 0.1 * nrow(dat)
            temp <- dat[c(order(dat[, 2], decreasing=TRUE)[1:nr],
                          order(dat[, 2], decreasing=FALSE)[1:nr]), ]
            y <- temp[, 2]
            x <- temp[, 1]
            fit <- stats::lm(y~x)
            graphics::abline(fit, lty=2, col="red")
            m <- summary(fit)$coefficients[2, 1]
            intercept <- summary(fit)$coefficients[1, 1]
            Rsq <- summary(fit)$r.squared
            posx <- graphics::par("usr")[1] + (graphics::par("usr")[2]
                                               - graphics::par("usr")[1]) * 0.15
            posy <- graphics::par("usr")[4] - (graphics::par("usr")[4]
                                               - graphics::par("usr")[3]) * 0.1
            labels <- base::paste("R? =", round(Rsq, digits=2), "\nslope =",
                                  round(m, digits=2))
            graphics::text(posx, posy, labels=labels)
            grDevices::dev.off()

            # Correction
            vals <- base::log2(10^feat_sample_mat_requant)
            features_sample_matrix_corrected <- base::as.data.frame(vals)

            # Correct slope
            correction <- 1 / m

            for(c in 1:ncol(features_sample_matrix_corrected)){
                val <- (features_sample_matrix_corrected[, c] * correction)
                features_sample_matrix_corrected[, c] <- val
            }

            return(list(features_sample_matrix_corrected, correction_data=dat,
                        correction_fit=fit, correction_factor=correction))
        }

        # Correction factor is supplied so correct data accordingly
        else{
            # Correction
            vals <- base::log2(10^feat_sample_mat_requant)
            features_sample_matrix_corrected <- base::as.data.frame(vals)

            for(c in 1:ncol(features_sample_matrix_corrected)){
                val <- features_sample_matrix_corrected[, c] * corr_factor
                features_sample_matrix_corrected[, c] <- val
            }

            return(list(features_sample_matrix_corrected,
                        correction_factor=corr_factor))
        }
    }

    Top3_Protein_Quant <- function(features, feat_sample_mat_requant,
                                   Alignment_scores=NULL, Quant_pvals=NULL,
                                   S2B=NULL, use_overlapping=TRUE, min_peps=2,
                                   quant_pvalue_cutoff=0.1,
                                   use_isotope_pmps=FALSE){
        label <- base::paste(round(0 / 1 * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title="Prepare Top3 quantification",
                                   label=label, min=0, max=1, width=300)
        close(pb)

        feat_sample_mat_requant <- base::as.data.frame(feat_sample_mat_requant)

        if(!is.null(Quant_pvals)){
            Quant_pvals <- base::as.data.frame(Quant_pvals)
        }

        if(!is.null(Alignment_scores)){
            Alignment_scores <- base::as.data.frame(Alignment_scores)
        }

        if(!is.null(S2B)){
            S2B <- base::as.data.frame(S2B)
        }

        # Top3 method
        # Input matrix with samples in cols and rows correspond to unique
        # peptides (log2 summed intensity over charge state and modification) of
        # a respective protein
        Top3_quant <- function(pep_matrix, features_temp,
                               Alignment_scores_temp=NULL,
                               Quant_pvals_temp=NULL, S2B_temp=NULL){
            sequence <- features_temp$Sequence

            if(length(sequence) > 0){
                pep_matrix <- stats::aggregate(2^pep_matrix,
                                               by=list(Sequence=sequence),
                                               FUN=sum, na.rm=TRUE)
                pep_matrix <- pep_matrix[, -1]
                pep_matrix[pep_matrix == 0] <- NA
                pep_matrix <- base::log2(pep_matrix)

                if(!is.null(Alignment_scores_temp)){
                    score_matrix <- stats::aggregate(Alignment_scores_temp,
                                                     by=list(Sequence=sequence),
                                                     FUN=mean, na.rm=TRUE)
                    score_matrix <- score_matrix[, -1]
                    score_matrix[score_matrix == 0] <- NA
                }

                else{
                    score_matrix <- matrix(nrow=nrow(pep_matrix),
                                           ncol=ncol(pep_matrix), NA)
                    colnames(score_matrix) <- colnames(pep_matrix)
                }

                if(!is.null(Quant_pvals_temp)){
                    pval_matrix <- stats::aggregate(Quant_pvals_temp,
                                                    by=list(Sequence=sequence),
                                                    FUN=mean, na.rm=TRUE)
                    pval_matrix <- pval_matrix[, -1]
                    pval_matrix[pval_matrix == 0] <- NA
                }

                else{
                    pval_matrix <- matrix(nrow=nrow(pep_matrix),
                                          ncol=ncol(pep_matrix), NA)
                    colnames(pval_matrix) <- colnames(pep_matrix)
                }

                if(!is.null(S2B_temp)){
                    S2B_matrix <- stats::aggregate(S2B_temp,
                                                   by=list(Sequence=sequence),
                                                   FUN=mean, na.rm=TRUE)
                    S2B_matrix <- S2B_matrix[, -1]
                    S2B_matrix[S2B_matrix == 0] <- NA
                }

                else{
                    S2B_matrix <- matrix(nrow=nrow(pep_matrix),
                                         ncol=ncol(pep_matrix), NA)
                    colnames(S2B_matrix) <- colnames(pep_matrix)
                }

                top3_res <- NULL
                top3_score_res <- NULL
                top3_pval_res <- NULL
                top3_S2B_res <- NULL

                for(c in 1:ncol(pep_matrix)){
                    top3 <- 2^pep_matrix[order(pep_matrix[, c], na.last=TRUE,
                                               decreasing=TRUE), c][1:3]
                    sums <- sum(top3, na.rm=TRUE)

                    if(!is.null(Alignment_scores_temp)){
                        ord <- order(pep_matrix[, c], na.last=TRUE,
                                     decreasing=TRUE)
                        top3_score <- score_matrix[ord, c][1:3]
                        median_score <- stats::median(top3_score, na.rm=TRUE)
                    }

                    else{
                        median_score <- NA
                    }

                    if(!is.null(Quant_pvals_temp)){
                        top3_pval <- pval_matrix[order(pep_matrix[, c],
                                                       na.last=TRUE,
                                                       decreasing=TRUE), c][1:3]
                        median_pval <- stats::median(top3_pval, na.rm=TRUE)
                    }

                    else{
                        median_pval <- NA
                    }

                    if(!is.null(S2B_temp)){
                        top3_S2B <- S2B_matrix[order(pep_matrix[, c],
                                                     na.last=TRUE,
                                                     decreasing=TRUE), c][1:3]
                        median_S2B <- stats::median(top3_S2B, na.rm=TRUE)
                    }

                    else{
                        median_S2B <- NA
                    }

                    if(sums == 0){
                        sums <- NA
                    }

                    else{
                        # Quantification only if at least 2 peptide
                        # quantifications are available
                        if(length(which(!is.na(top3))) >= min_peps){
                            sums <- sums / length(which(!is.na(top3)))
                        }

                        else{
                            sums <- NA
                            median_score <- NA
                            median_pval <- NA
                            median_S2B <- NA
                        }
                    }

                    top3_res <- append(top3_res, base::log2(sums))
                    top3_score_res <- append(top3_score_res, median_score)
                    top3_pval_res <- append(top3_pval_res, median_pval)
                    top3_S2B_res <- append(top3_S2B_res, median_S2B)
                }
            }

            else{
                top3_res <- rep(NA, ncol(pep_matrix))
                top3_score_res <- rep(NA, ncol(score_matrix))
                top3_pval_res <- rep(NA, ncol(pval_matrix))
                top3_S2B_res <- rep(NA, ncol(pval_matrix))
            }

            names(top3_res) <- colnames(pep_matrix)
            names(top3_score_res) <- base::paste(colnames(score_matrix),
                                                 "_median_score", sep="")
            names(top3_pval_res) <- base::paste(colnames(pval_matrix),
                                                "_median_pvals", sep="")
            names(top3_S2B_res) <- base::paste(colnames(pval_matrix),
                                               "_median_S2B", sep="")

            return(append(top3_res,
                          append(top3_score_res,
                                 append(top3_pval_res, top3_S2B_res))))
        }

        c1 <- features$Protein != ""
        c2 <- !grepl(";", features$Sequence)
        c3 <- !grepl("_i|_pmp", features$Feature_name)

        if(use_isotope_pmps == FALSE){
            loc <- which(c1 & c2 & c3)
            features_temp <- features[loc, ]
            feat_sample_mat_requant_temp <- feat_sample_mat_requant[loc, ]
            Quant_pvals_temp <- Quant_pvals[loc, ]
            Alignment_scores_temp <- Alignment_scores[loc, ]
            S2B_temp <- S2B[loc, ]
        }

        else{
            loc <- which(c1 & c2)
            features_temp <- features[loc, ]
            feat_sample_mat_requant_temp <- feat_sample_mat_requant[loc, ]
            Quant_pvals_temp <- Quant_pvals[loc, ]
            Alignment_scores_temp <- Alignment_scores[loc, ]
            S2B_temp <- S2B[loc, ]
        }

        if(!is.na(quant_pvalue_cutoff)){
            rmins <- matrixStats::rowMins(as.matrix(Quant_pvals_temp),
                                          na.rm=TRUE)
            selection <- which(rmins < quant_pvalue_cutoff)

            features_temp <- features[selection, ]
            feat_sample_mat_requant_temp <- feat_sample_mat_requant[selection, ]
            Quant_pvals_temp <- Quant_pvals[selection, ]
            Alignment_scores_temp <- Alignment_scores[selection, ]
            S2B_temp <- S2B[selection, ]
        }

        if(nrow(features_temp) > 0){
            # Use also features where a peptide is shared between 2 or more
            # proteins
            if(use_overlapping == TRUE){
                prts <- as.character(stringr::str_split(features_temp$Protein,
                                                        "\\||;", simplify=TRUE))
                unique_proteins <- sort(unique(prts))
            }

            else{
                loc <- which(!grepl("\\||;", features_temp$Protein))
                prts <- features_temp$Protein[loc]
                unique_proteins <- sort(unique(prts))
            }

            if(any(unique_proteins == "")){
                loc <- which(unique_proteins == "")
                unique_proteins <- unique_proteins[-loc]
            }

            mat <- matrix(ncol=1 + (4 * ncol(feat_sample_mat_requant)),
                          nrow=length(unique_proteins), 0)
            protein_TOP3 <- base::as.data.frame(mat)
            rownames(protein_TOP3) <- unique_proteins
            mas <- base::paste("median_alignment_score_",
                               colnames(feat_sample_mat_requant), sep="")
            mqp <- base::paste("median_quant_pvals_",
                               colnames(feat_sample_mat_requant), sep="")
            ms2b <- base::paste("median_S2B_",
                                colnames(feat_sample_mat_requant), sep="")
            colnames(protein_TOP3) <- c("num_quant_features",
                                        colnames(feat_sample_mat_requant),
                                        mas, mqp, ms2b)

            max <- nrow(protein_TOP3)
            label <- base::paste(round(0 / max * 100, 0), "% done")
            pb <- tcltk::tkProgressBar(title="Perform Top3 quantification",
                                       label=label, min=0, max=max, width=300)
            start_time <- Sys.time()
            updatecounter <- 0
            time_require <- 0

            for(i in 1:length(unique_proteins)){
                if(use_overlapping == TRUE){
                    ind <- which(grepl(unique_proteins[i],
                                       features_temp$Protein))
                }

                else{
                    ind <- which(features_temp$Protein == unique_proteins[i])
                }

                if(length(ind) > 0){
                    res <- Top3_quant(
                        pep_matrix=feat_sample_mat_requant_temp[ind, ],
                        features_temp=features_temp[ind, ],
                        Alignment_scores_temp=Alignment_scores[ind, ],
                        Quant_pvals_temp=Quant_pvals[ind, ],
                        S2B_temp=S2B_temp[ind, ]
                    )
                    data.table::set(protein_TOP3, as.integer(i),
                                    as.integer(1:ncol(protein_TOP3)),
                                    value=as.list(c(length(ind),
                                                    as.numeric(res))))
                }

                else{
                    data.table::set(protein_TOP3, as.integer(i), as.integer(1),
                                    value=0)
                }

                updatecounter <- updatecounter + 1

                if(updatecounter >= 10){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                                  digits=0))

                    updatecounter <- 0
                    label <- base::paste(round(i / max * 100, 0), " % done (",
                                         i, "/", max, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }

            close(pb)
        }

        else{
            ncols <- 1 + (4 * ncol(feat_sample_mat_requant))
            protein_TOP3 <- base::as.data.frame(matrix(ncol=ncols, nrow=0, 0))
            cnames <- c("num_quant_features", colnames(feat_sample_mat_requant),
                        base::paste("median_alignment_score_",
                                    colnames(feat_sample_mat_requant), sep=""),
                        base::paste("median_quant_pvals_",
                                    colnames(feat_sample_mat_requant), sep=""),
                        base::paste("median_S2B_",
                                    colnames(feat_sample_mat_requant), sep=""))
            colnames(protein_TOP3) <- cnames
        }

        return(protein_TOP3)
    }

    Total_Protein_Quant <- function(features, feat_sample_mat_requant,
                                    Alignment_scores=NULL, Quant_pvals=NULL,
                                    S2B=NULL, use_overlapping=TRUE, min_peps=2,
                                    quant_pvalue_cutoff=0.1,
                                    use_isotope_pmps=FALSE){
        label <- base::paste(round(0 / 1 * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title="Prepare Total quantification",
                                   label=label, min=0, max=1, width=300)
        close(pb)

        feat_sample_mat_requant <- base::as.data.frame(feat_sample_mat_requant)

        if(!is.null(Quant_pvals)){
            Quant_pvals <- base::as.data.frame(Quant_pvals)
        }

        if(!is.null(Alignment_scores)){
            Alignment_scores <- base::as.data.frame(Alignment_scores)
        }

        if(!is.null(S2B)){
            S2B <- base::as.data.frame(S2B)
        }

        # Total method
        # Input matrix with samples in cols and rows correspond to unique
        # peptides (log2 summed intensity over charge state and modification) of
        # a respective protein
        Total_quant <- function(pep_matrix, features_temp,
                                Alignment_scores_temp=NULL,
                                Quant_pvals_temp=NULL, S2B_temp=NULL,
                                Quant_cutoff=4){
            sequence <- features_temp$Sequence

            if(length(sequence) > 0){
                pep_matrix <- stats::aggregate(2^pep_matrix,
                                               by=list(Sequence=sequence),
                                               FUN=sum, na.rm=TRUE)
                pep_matrix <- pep_matrix[, -1]
                pep_matrix[pep_matrix == 0] <- NA
                pep_matrix <- base::log2(pep_matrix)

                if(!is.null(Alignment_scores_temp)){
                    score_matrix <- stats::aggregate(Alignment_scores_temp,
                                                     by=list(Sequence=sequence),
                                                     FUN=mean, na.rm=TRUE)
                    score_matrix <- score_matrix[, -1]
                    score_matrix[score_matrix == 0] <- NA
                }

                else{
                    score_matrix <- matrix(nrow=nrow(pep_matrix),
                                           ncol=ncol(pep_matrix), NA)
                    colnames(score_matrix) <- colnames(pep_matrix)
                }

                if(!is.null(Quant_pvals_temp)){
                    pval_matrix <- stats::aggregate(Quant_pvals_temp,
                                                    by=list(Sequence=sequence),
                                                    FUN=mean, na.rm=TRUE)
                    pval_matrix <- pval_matrix[, -1]
                    pval_matrix[pval_matrix == 0] <- NA
                }

                else{
                    pval_matrix <- matrix(nrow=nrow(pep_matrix),
                                          ncol=ncol(pep_matrix), NA)
                    colnames(pval_matrix) <- colnames(pep_matrix)
                }

                if(!is.null(S2B_temp)){
                    S2B_matrix <- stats::aggregate(S2B_temp,
                                                   by=list(Sequence=sequence),
                                                   FUN=mean, na.rm=TRUE)
                    S2B_matrix <- S2B_matrix[, -1]
                    S2B_matrix[S2B_matrix == 0] <- NA
                }

                else{
                    S2B_matrix <- matrix(nrow=nrow(pep_matrix),
                                         ncol=ncol(pep_matrix), NA)
                    colnames(S2B_matrix) <- colnames(pep_matrix)
                }

                total_res <- NULL
                total_score_res <- NULL
                total_pval_res <- NULL
                total_S2B_res <- NULL

                for(c in 1:ncol(pep_matrix)){
                    total <- 2^pep_matrix[order(pep_matrix[, c], na.last=TRUE,
                                                decreasing=TRUE), c]
                    sums <- sum(total, na.rm=TRUE)

                    if(!is.null(Alignment_scores_temp)){
                        total_score <- score_matrix[order(pep_matrix[, c],
                                                          na.last=TRUE,
                                                          decreasing=TRUE), c]
                        median_score <- stats::weighted.mean(total_score,
                                                             total, na.rm=TRUE)
                    }

                    else{
                        median_score <- NA
                    }

                    if(!is.null(Quant_pvals_temp)){
                        total_pval <- pval_matrix[order(pep_matrix[, c],
                                                        na.last=TRUE,
                                                        decreasing=TRUE), c]
                        median_pval <- stats::weighted.mean(total_pval,total,
                                                            na.rm=TRUE)
                    }

                    else{
                        median_pval <- NA
                    }

                    if(!is.null(S2B_temp)){
                        total_S2B <- S2B_matrix[order(pep_matrix[, c],
                                                      na.last=TRUE,
                                                      decreasing=TRUE), c]
                        median_S2B <- stats::weighted.mean(total_S2B, total,
                                                           na.rm=TRUE)
                    }

                    else{
                        median_S2B <- NA
                    }

                    if(sums == 0){
                        sums <- NA
                    }

                    else{
                        # Quantification only if at least n peptide
                        # quantifications are available
                        if(length(which(!is.na(total))) >= min_peps){
                            sums <- sums / length(which(!is.na(total)))
                        }

                        else{
                            sums <- NA
                            median_score <- NA
                            median_pval <- NA
                            median_S2B <- NA
                        }
                    }

                    total_res <- append(total_res, base::log2(sums))
                    total_score_res <- append(total_score_res, median_score)
                    total_pval_res <- append(total_pval_res, median_pval)
                    total_S2B_res <- append(total_S2B_res, median_S2B)
                }
            }

            else{
                total_res <- rep(NA, ncol(pep_matrix))
                total_score_res <- rep(NA, ncol(score_matrix))
                total_pval_res <- rep(NA, ncol(pval_matrix))
                total_S2B_res <- rep(NA, ncol(pval_matrix))
            }

            names(total_res) <- colnames(pep_matrix)
            names(total_score_res) <- base::paste(colnames(score_matrix),
                                                  "_median_score", sep="")
            names(total_pval_res) <- base::paste(colnames(pval_matrix),
                                                 "_median_pvals", sep="")
            names(total_S2B_res) <- base::paste(colnames(pval_matrix),
                                                "_median_S2B", sep="")

            return(append(total_res,
                          append(total_score_res,
                                 append(total_pval_res, total_S2B_res))))
        }

        c1 <- features$Protein != ""
        c2 <- !grepl(";", features$Sequence)
        c3 <- !grepl("_i|_pmp", features$Feature_name)

        if(use_isotope_pmps == FALSE){
            loc <- which(c1 & c2 & c3)

            features_temp <- features[loc, ]
            feat_sample_mat_requant_temp <- feat_sample_mat_requant[loc, ]
            Quant_pvals_temp <- Quant_pvals[loc, ]
            Alignment_scores_temp <- Alignment_scores[loc, ]
            S2B_temp <- S2B[loc, ]
        }

        else{
            loc <- which(c1 & c2)

            features_temp <- features[loc, ]
            feat_sample_mat_requant_temp <- feat_sample_mat_requant[loc, ]
            Quant_pvals_temp <- Quant_pvals[loc, ]
            Alignment_scores_temp <- Alignment_scores[loc, ]
            S2B_temp <- S2B[loc, ]
        }

        if(!is.na(quant_pvalue_cutoff)){
            mns <- matrixStats::rowMins(as.matrix(Quant_pvals_temp), na.rm=TRUE)
            selection <- which(mns < quant_pvalue_cutoff)

            features_temp <- features[selection, ]
            feat_sample_mat_requant_temp <- feat_sample_mat_requant[selection, ]
            Quant_pvals_temp <- Quant_pvals[selection, ]
            Alignment_scores_temp <- Alignment_scores[selection, ]
            S2B_temp <- S2B[selection, ]
        }

        if(nrow(features_temp) > 0){
            # Use also features where a peptide is shared between 2 or more
            # proteins
            if(use_overlapping == TRUE){
                ssplit <- stringr::str_split(features_temp$Protein, "\\||;",
                                             simplify=TRUE)
                unique_proteins <- sort(unique(as.character(ssplit)))
            }

            else{
                loc <- which(!grepl("\\||;", features_temp$Protein))
                unique_proteins <- sort(unique(features_temp$Protein[loc]))
            }

            if(any(unique_proteins == "")){
                loc <- which(unique_proteins == "")
                unique_proteins <- unique_proteins[-loc]
            }

            ncols <- 1 + (4 * ncol(feat_sample_mat_requant))
            mat <- matrix(ncol=ncols, nrow=length(unique_proteins), 0)
            protein_total <- base::as.data.frame(mat)
            rownames(protein_total) <- unique_proteins
            cnames <- c("num_quant_features", colnames(feat_sample_mat_requant),
                        base::paste("median_alignment_score_",
                                    colnames(feat_sample_mat_requant), sep=""),
                        base::paste("median_quant_pvals_",
                                    colnames(feat_sample_mat_requant), sep=""),
                        base::paste("median_S2B_",
                                    colnames(feat_sample_mat_requant), sep=""))
            colnames(protein_total) <- cnames

            max <- nrow(protein_total)
            label <- base::paste(round(0 / max * 100, 0), "% done")
            pb <- tcltk::tkProgressBar(title="Perform total quantification",
                                       label=label, min=0, max=max, width=300)
            start_time <- Sys.time()
            updatecounter <- 0
            time_require <- 0
            for(i in 1:length(unique_proteins)){
                if(use_overlapping == TRUE){
                    ind <- which(grepl(unique_proteins[i],
                                       features_temp$Protein))
                }

                else{
                    ind <- which(features_temp$Protein == unique_proteins[i])
                }

                if(length(ind) > 0){
                    res <- Total_quant(
                        pep_matrix=feat_sample_mat_requant_temp[ind, ],
                        features_temp=features_temp[ind, ],
                        Alignment_scores_temp=Alignment_scores[ind, ],
                        Quant_pvals_temp=Quant_pvals[ind, ],
                        S2B_temp=S2B_temp[ind, ]
                    )
                    data.table::set(protein_total, as.integer(i),
                                    as.integer(1:ncol(protein_total)),
                                    value=as.list(c(length(ind),
                                                    as.numeric(res))))
                }

                else{
                    data.table::set(protein_total, as.integer(i), as.integer(1),
                                    value=0)
                }

                updatecounter <- updatecounter + 1
                if(updatecounter >= 10){
                    time_elapsed <- difftime(Sys.time(), start_time,
                                             units="secs")
                    time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
                    td <- lubridate::seconds_to_period(time_require)
                    time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                            lubridate::minute(td),
                                            round(lubridate::second(td),
                                            digits=0))

                    updatecounter <- 0
                    label <- base::paste(round(i / max * 100, 0), " % done (",
                                         i, "/", max, ", Time require: ",
                                         time_require, ")", sep="")
                    tcltk::setTkProgressBar(pb, i, label=label)
                }
            }

            close(pb)
        }

        else{
            ncols <- 1 + (4 * ncol(feat_sample_mat_requant))
            protein_total <- base::as.data.frame(matrix(ncol=ncols, nrow=0, 0))
            cnames <- c("num_quant_features", colnames(feat_sample_mat_requant),
                        base::paste("median_alignment_score_",
                                    colnames(feat_sample_mat_requant), sep=""),
                        base::paste("median_quant_pvals_",
                                    colnames(feat_sample_mat_requant), sep=""),
                        base::paste("median_S2B_",
                                    colnames(feat_sample_mat_requant), sep=""))
            colnames(protein_total) <- cnames
        }

        return(protein_total)
    }

    crap <- gc(FALSE)

    # Perform protein level aggregation
    if(!is.null(path_to_MaxQ_output)){
        print(paste0(Sys.time()," Perform protein-level aggregation"))

        # Load MaxQ peptide results
        pth <- base::paste(path_to_MaxQ_output, "/peptides.txt", sep="")
        MaxQ_peptides <- utils::read.table(pth, sep="\t", header=TRUE)
        loc1 <- which(MaxQ_peptides$Potential.contaminant == "+"
                     | MaxQ_peptides$Reverse == "+")
        MaxQ_peptides <- MaxQ_peptides[-loc1,]
        loc2 <- which(grepl("Intensity\\.", colnames(MaxQ_peptides)))
        MaxQ_peptides_quant <- MaxQ_peptides[, loc2]
        MaxQ_peptides_quant[MaxQ_peptides_quant == 0] <- NA
        MaxQ_peptides_quant <- base::log2(MaxQ_peptides_quant)
        df <- base::data.frame(
            Sequence=MaxQ_peptides$Sequence,
            Leading_razor=MaxQ_peptides$Leading.razor.protein
        )
        MaxQ_peptides_leading_razor <- df
        lrzr <- as.character(MaxQ_peptides_leading_razor$Leading_razor)
        MaxQ_peptides_leading_razor$Leading_razor <- lrzr
        MaxQ_peptides <- base::data.frame(Sequence=MaxQ_peptides$Sequence,
                                          MaxQ_peptides_quant)
        pth <- base::paste(path_to_MaxQ_output, "/proteinGroups.txt", sep="")
        MaxQ_protein_groups <- utils::read.table(pth, sep="\t", header=TRUE)

        # Check if IDs were correctly parsed, if not, try to parse with
        # SwissProt or Trembl
        # Not correctly parsed but contains trembl or swissprot fasta headers
        c1 <- !any(colnames(MaxQ_protein_groups) == "Gene.names")
        c2 <- any(colnames(MaxQ_protein_groups) == "Protein.IDs")
        if(c1 & c2){
            print(paste0(Sys.time(), " Fasta file was not correctly parsed ",
                         "during search. Try to base::paste fasta headers ..."))
            print(paste0(Sys.time(), " Detected Swiss-Prot and/or TrEMBL ",
                         "fasta headers."))
            if(any(grepl(">sp|>tr", MaxQ_protein_groups$Fasta.headers))){
                # Protein level
                # Parsing was not performed correctly so we have to try to do
                # this here expecting swissprot or trembl fasta headers
                MaxQ_protein_groups$Gene.names <- ""
                MaxQ_protein_groups$Organism <- ""
                pids <- as.character(MaxQ_protein_groups$Protein.IDs)
                MaxQ_protein_groups$Protein.IDs <- pids

                hdr <- base::gsub(">", "", MaxQ_protein_groups$Fasta.headers)
                fasta_headers <- stringr::str_split(hdr, "\\|", simplify=TRUE)

                # Extract gene name and species information
                GN <- vector("character", nrow(MaxQ_protein_groups))
                ID <- vector("character", nrow(MaxQ_protein_groups))
                regx <- gregexpr("GN=", MaxQ_protein_groups$Fasta.headers)
                GN_start <- unlist(lapply(base::paste(regx, sep=","), `[[`, 1))
                GN_start <- base::gsub("c\\(|\\)", "", GN_start)

                for(i in 1:nrow(MaxQ_protein_groups)){
                    if(GN_start[i] != "-1"){
                        ssplit <- stringr::str_split(GN_start[i], ",")
                        indices <- as.numeric(unlist(ssplit)) + 3
                        sstr <- substring(MaxQ_protein_groups$Fasta.headers[i],
                                          indices, indices + 10)
                        stop <- regexpr(" ", sstr)
                        GN_temp <- substring(
                            MaxQ_protein_groups$Fasta.headers[i],
                            indices,
                            indices + stop - 2
                        )
                        GN[i] <- base::paste(GN_temp, collapse=";")
                    }

                    loc <- which(fasta_headers[i, c(2, 4, 6)] != "")
                    header <- fasta_headers[i, c(2, 4, 6)][loc]

                    if(length(header) > 0){
                        ID[i] <- base::paste(header, collapse=";")
                    }
                }

                pids <- MaxQ_protein_groups$Protein.IDs
                CON_REV <- !grepl("^CON_|^REV_",
                                  MaxQ_protein_groups$Protein.IDs)
                npids <- ifelse(CON_REV, ID, pids)
                MaxQ_protein_groups$Protein.IDs <- npids
                MaxQ_protein_groups$Majority.protein.IDs <- npids
                MaxQ_protein_groups$Gene.names <- ifelse(CON_REV, GN, "")

                # Peptide level
                cnd <- grepl("\\|", MaxQ_peptides_leading_razor$Leading_razor)

                if(any(cnd)){
                    lrzr <- MaxQ_peptides_leading_razor$Leading_razor
                    temp <- stringr::str_split(lrzr, "\\|", simplify=TRUE)
                    res <- ifelse(temp[, 2] != "", temp[, 2], lrzr)
                    MaxQ_peptides_leading_razor$Leading_razor <- res
                }

                # Requantification features
                temp <- stringr::str_split(features$Protein, "\\||;",
                                           simplify=TRUE)
                ID <- vector("character", nrow(temp))

                for(i in 1:nrow(temp)){
                    sel <- which(temp[i, ] == "sp") + 1

                    if(length(sel) > 0){
                        ID[i] <- base::paste(temp[i, sel], collapse=";")
                    }

                    else{
                        sel <- which(temp[i, ] != "")
                        ID[i] <- base::paste(temp[i, sel], collapse=";")
                    }
                }

                features$Protein <- ID
            }

            else{
                print(paste0(Sys.time(), " No supported fasta headers were ",
                             "detected. Next steps might be not fully ",
                             "working."))
            }
        }

        # Correct for overestimation of high intense features
        # No correction
        if(abundance_estimation_correction == FALSE){
            feat_w_backgr_int <- base::log2(10^feat_w_backgr_int)
            aux <- 10^feat_w_backgr_int_imputed
            feat_w_backgr_int_imputed <- base::log2(aux)
        }

        # Correct abundance estimations based on MaxQ peptide abundance
        # estimations
        else{
            # Determine abundance correction factors based on MaxQ peptide
            # intensities and data without imputation
            cor_res <- correct_intensities(
                features,
                feat_sample_mat_requant=feat_w_backgr_int,
                pval_quant=pval_sig_w_bckgrnd_quant,
                MaxQ_peptides_quant=MaxQ_peptides,
                main="Signal_Background_intensity"
            )

            feat_w_backgr_int <- cor_res[[1]]
            lst <- list(correction_data=cor_res[[2]],
                        correction_fit=cor_res[[3]],
                        correction_factor=cor_res[[4]])
            QC_data[["Abundance_correction"]] <- lst

            cor_res <- correct_intensities(
                features,
                feat_sample_mat_requant=feat_w_backgr_int_imputed,
                pval_quant=pval_sig_w_bckgrnd_quant,
                MaxQ_peptides_quant=MaxQ_peptides,
                main="Signal_Background_intensity_imputed",
                corr_factor=cor_res[[4]]
            )
            feat_w_backgr_int_imputed <- cor_res[[1]]
        }

        # Get leading razor ID per peptide sequence
        features$Protein <- as.character(features$Protein)
        features$all_matching_Proteins <- features$Protein
        loc <- match(features$Sequence, MaxQ_peptides_leading_razor$Sequence)
        aux1 <- MaxQ_peptides_leading_razor$Leading_razor[loc]
        features$Protein <- as.character(aux1)
        aux2 <- features$all_matching_Proteins[is.na(features$Protein)]
        features$Protein[is.na(features$Protein)] <- aux2

        # Perform peptide level LFQ
        if(calc_peptide_LFQ == TRUE){
            LFQ_peptide_quant_process <- function(
                features,
                features_quant,
                n_cores,
                label="Perform peptide-LFQ-quantification",
                num_ratio_samples=NA,
                TopN=5,
                seed=1
            ){
                set.seed(seed)
                # Prepare for MaxLFQ algorithm
                calculate_LFQ <- function(peptide_quant_data, min_num_ratios=2,
                                          num_ratio_samples=NA){

                    # Error function which is used to optimize ratios
                    least_square_error <- function(par, ratio_mat){
                        sum <- 0
                        vals <- NULL

                        if(ncol(ratio_mat) > 1){
                            for(c in 1:(ncol(ratio_mat) - 1)){
                                for(r in (c + 1):nrow(ratio_mat)){
                                    val <- (base::log2(ratio_mat[r, c])
                                            - base::log2(par[r])
                                            + base::log2(par[c]))^2
                                    vals <- append(vals, val)

                                    if(!is.na(val)){
                                        sum <- sum + val
                                    }
                                }
                            }
                        }

                        if(sum == 0){
                            sum <- NA
                        }

                        return(sum)
                    }

                    # Unlog intensities
                    peptide_quant_data <- 2^peptide_quant_data
                    # Calculate summed intensities per sample
                    totalsum_per_sample <- colSums(peptide_quant_data,
                                                   na.rm=TRUE)
                    # Determine median ratio matrix between all samples
                    mat <- matrix(nrow=ncol(peptide_quant_data),
                                  ncol=ncol(peptide_quant_data))
                    ratio_mat <- base::as.data.frame(mat)

                    if(ncol(ratio_mat) > 1){
                        for(c in 1:(ncol(ratio_mat) - 1)){
                            for(r in (c + 1):nrow(ratio_mat)){
                                ratios <- (peptide_quant_data[, r]
                                           / peptide_quant_data[, c])

                                cnd <- length(which(!is.na(ratios)))

                                if(cnd >= min_num_ratios){
                                    ratio_mat[r, c] <- stats::median(ratios,
                                                                     na.rm=TRUE)
                                }
                            }
                        }
                    }

                    # Calculate ratios over all samples
                    if(is.na(num_ratio_samples)){
                        # Define start parameter
                        start_par <- c(rep(1, ncol(ratio_mat)))
                        # Now find optimum in ratios to best recover true
                        # observed ratios between samples
                        res_ratio <- NULL
                        try(res_ratio <- stats::optim(par=start_par,
                                                      fn=least_square_error,
                                                      ratio_mat=ratio_mat,
                                                      lower=1, upper=100,
                                                      method="L-BFGS-B"),
                            silent=TRUE)

                        if(!is.null(res_ratio)){
                            # Normalize ratios to sample with highest intensity
                            mx <- max(totalsum_per_sample)
                            loc <- which(totalsum_per_sample == mx)
                            ratio_norm <- res_ratio$par / res_ratio$par[loc]
                            # Finally calculate log2 lfq protein intensities per
                            # sample
                            lfq <- base::log2(ratio_norm
                                              * totalsum_per_sample[loc])
                            # Remove quant values for samples were no ratios
                            # were available
                            c1 <- colSums(!is.na(ratio_mat)) == 0
                            c2 <- rowSums(!is.na(ratio_mat)) == 0

                            if(any(c1 & c2)){
                                sel <- as.numeric(which(c1 & c2))
                                lfq[sel] <- NA
                            }
                        }

                        else{
                            lfq <- rep(NA, ncol(peptide_quant_data))
                        }
                    }

                    # Calculate ratios only for a subset of samples
                    else{
                        mat <- matrix(nrow=ncol(peptide_quant_data),
                                      ncol=ncol(peptide_quant_data), 0)
                        lfq <- base::as.data.frame(mat)

                        for(s in 1:ncol(peptide_quant_data)){
                            # Randomly select up to 6 other samples from list of
                            # samples which also contain quantifications
                            csum <- colSums(peptide_quant_data, na.rm=TRUE)
                            samples_with_quant <- as.numeric(which(csum > 0))
                            loc <- which(samples_with_quant != s)
                            samples_with_quant <- samples_with_quant[loc]

                            if(length(samples_with_quant) > 0){
                                len <- length(samples_with_quant)
                                samps_comparison_x_samp <- sample(
                                    samples_with_quant,
                                    ifelse(len > num_ratio_samples,
                                           num_ratio_samples, len)
                                )
                            }

                            else{
                                samps_comparison_x_samp <- sample(
                                    c(1:ncol(peptide_quant_data))[-i],
                                    size=num_ratio_samples,
                                    replace=FALSE
                                )
                            }

                            loc <- c(s, sort(samps_comparison_x_samp))
                            ratio_mat_temp <- ratio_mat[loc, loc]

                            # To few of selected random samples show an observed
                            # intensity ratio but protein is quantified in
                            # current sample
                            c1 <- length(which(!is.na(ratio_mat_temp))) < 3
                            c2 <- any(!is.na(peptide_quant_data[, s]))

                            # XXX: Why empty?
                            if(c1 & c2){
                                # Randomly select other samples
                            }

                            totalsum_per_sample_temp <- totalsum_per_sample[loc]
                            # Define start parameter
                            start_par <- c(rep(1, ncol(ratio_mat_temp)))
                            # Now find optimum in ratios to best recover true
                            # observed ratios between samples
                            res_ratio <- NULL
                            try(res_ratio <- stats::optim(
                                par=start_par,
                                fn=least_square_error,
                                ratio_mat=ratio_mat_temp,
                                lower=1,
                                upper=100,
                                method="L-BFGS-B"
                            ), silent=TRUE)

                            if(!is.null(res_ratio)){
                                # Normalize ratios to sample with highest
                                # intensity
                                mx <- max(totalsum_per_sample_temp, na.rm=TRUE)
                                loc <- which(totalsum_per_sample_temp == mx)
                                ratio_norm <- res_ratio$par / res_ratio$par[loc]
                                # Finally calculate log2 lfq protein intensities
                                # per sample
                                ts <- totalsum_per_sample_temp[loc]
                                prod <- ratio_norm * ts
                                temp_lfq <- base::log2(prod)
                                lst <- c(s, sort(samps_comparison_x_samp))
                                val <- as.integer(lst)
                                data.table::set(lfq, as.integer(s), val,
                                                as.list(temp_lfq))
                            }
                        }

                        lfq[lfq == 0] <- NA
                        lfq <- matrixStats::colMedians(as.matrix(lfq),
                                                       na.rm=TRUE)

                        # Remove quant values for samples were no ratios were
                        # available
                        c1 <- colSums(!is.na(ratio_mat)) == 0
                        c2 <- rowSums(!is.na(ratio_mat)) == 0

                        if(any(c1 & c2)){
                            sel <- as.numeric(which(c1 & c2))
                            lfq[sel] <- NA
                        }
                    }

                    return(lfq)
                }

                # txtProgressBar from package pbarETA (Francesco Napolitano)
                # License: LGPL-3
                txtProgressBar <- function(min=0, max=1, initial=0, char="=",
                                           width=NA, title, label, style=3,
                                           file=""){
                    formatTime <- function(seconds){
                        if(seconds == Inf || is.nan(seconds) || is.na(seconds)){
                            return("NA")
                        }

                        seconds <- round(seconds)
                        sXmin <- 60
                        sXhr <- sXmin * 60
                        sXday <- sXhr * 24
                        sXweek <- sXday * 7
                        sXmonth <- sXweek * 4.22
                        sXyear <- sXmonth * 12
                        years <- floor(seconds / sXyear)
                        seconds <- seconds - years * sXyear
                        months <- floor(seconds / sXmonth)
                        seconds <- seconds - months * sXmonth
                        weeks <- floor(seconds / sXweek)
                        seconds <- seconds - weeks * sXweek
                        days <- floor(seconds / sXday)
                        seconds <- seconds - days * sXday
                        hours <- floor(seconds / sXhr)
                        seconds <- seconds - hours * sXhr
                        minutes <- floor(seconds / sXmin)
                        seconds <- seconds - minutes * sXmin
                        ETA <- c(years, months, days, hours, minutes, seconds)
                        startst <- which(ETA > 0)[1]

                        if(is.na(startst)){
                            startst <- 6
                        }

                        starts <- min(startst, 4)
                        len <- length(ETA)
                        fmtstr <- rep("%02d", length(ETA))[startst:len]
                        fmtstr <- base::paste(fmtstr, collapse=":")
                        lst <- as.list(c(as.list(fmtstr), ETA[startst:len]))

                        return(do.call(sprintf, lst))
                    }

                    c1 <- !identical(file, "")
                    c2 <- !(inherits(file, "connection") && isOpen(file))

                    if(c1 && c2){
                        stop("'file' must be \"\" or an open connection object")
                    }

                    if(!style %in% 1L:3L){
                        style <- 1
                    }

                    .val <- initial
                    .killed <- FALSE
                    .nb <- 0L
                    .pc <- -1L
                    .time0 <- NA
                    .timenow <- NA
                    .firstUpdate <- T
                    nw <- nchar(char, "w")

                    if(is.na(width)){
                        width <- getOption("width")

                        if(style == 3L){
                            width <- width - 10L
                        }

                        width <- trunc(width / nw)
                    }

                    if(max <= min){
                        stop("must have 'max' > 'min'")
                    }

                    up1 <- function(value){
                        if(!is.finite(value) || value < min || value > max){
                            return()
                        }

                        .val <<- value
                        nb <- round(width * (value - min) / (max - min))

                        if(.nb < nb){
                            cat(base::paste(rep.int(char, nb - .nb),
                                            collapse=""), file=file)
                            utils::flush.console()
                        }

                        else if(.nb > nb){
                            cat("\r", base::paste(rep.int(" ", .nb * nw),
                                                  collapse=""), "\r",
                                base::paste(rep.int(char, nb), collapse=""),
                                sep="", file=file)
                            utils::flush.console()
                        }

                        .nb <<- nb
                    }

                    up2 <- function(value){
                        if(!is.finite(value) || value < min || value > max){
                            return()
                        }

                        .val <<- value
                        nb <- round(width * (value - min) / (max - min))

                        if(.nb <= nb){
                            cat("\r", base::paste(rep.int(char, nb),
                                                  collapse=""), sep="",
                                file=file)
                            utils::flush.console()
                        }

                        else{
                            cat("\r", base::paste(rep.int(" ", .nb * nw),
                                                  collapse=""), "\r",
                                base::paste(rep.int(char, nb), collapse=""),
                                            sep="", file=file)
                            utils::flush.console()
                        }

                        .nb <<- nb
                    }

                    up3 <- function(value, calledOnCreation=FALSE){
                        timenow <- proc.time()[["elapsed"]]

                        if(!calledOnCreation && .firstUpdate){
                            .time0 <<- timenow
                            .timenow <<- timenow
                            .firstUpdate <<- FALSE
                        }

                        if(!is.finite(value) || value < min || value > max){
                            return()
                        }

                        .val <<- value
                        nb <- round(width * (value - min) / (max - min))
                        pc <- round(100 * (value - min) / (max - min))

                        if (nb == .nb && pc == .pc && timenow - .timenow < 1){
                            return()
                        }

                        .timenow <<- timenow
                        span <- timenow - .time0
                        timeXiter <- span / (.val - min)
                        ETA <- (max - .val) * timeXiter
                        ETAstr <- formatTime(ETA)
                        ri <- rep.int(" ", nw * width + 6)
                        cat(base::paste(c("\r  |", ri), collapse=""), file=file)
                        sf <- sprintf("| %3d%%", pc), ", ETA ", ETAstr)
                        cat(base::paste(c("\r  |", rep.int(char, nb),
                                          rep.int(" ", nw * (width - nb)), sf,
                                        collapse=""), file=file)
                        utils::flush.console()
                        .nb <<- nb
                        .pc <<- pc
                    }

                    getVal <- function(){
                        .val
                    }

                    kill <- function(){
                        if (!.killed) {
                            cat("\n", file=file)
                            utils::flush.console()
                            .killed <<- TRUE
                        }
                    }

                    up <- switch(style, up1, up2, up3)
                    up(initial, TRUE)
                    structure(list(getVal=getVal, up=up, kill=kill),
                              class="txtProgressBar")
                }

                sel <- which(!grepl(";|,", features$Sequence))
                aux <- base::paste(features$Sequence[sel],
                                   features$Modifications[sel], sep="_")
                unique_peptides <- sort(unique(aux))
                rgx <- regexpr("_", unique_peptides)
                unique_seq <- base::substr(unique_peptides, 1, rgx - 1)
                unique_mod <- base::substr(unique_peptides, rgx + 1,
                                           nchar(unique_peptides))

                mat <- matrix(ncol=ncol(features_quant),
                              nrow=length(unique_peptides), 0)
                LFQ_peptide_quant <- base::as.data.frame(mat)
                colnames(LFQ_peptide_quant) <- colnames(features_quant)

                # Perform LFQ quantification
                cl <- snow::makeCluster(n_cores)
                doSNOW::registerDoSNOW(cl)
                iterations <- nrow(LFQ_peptide_quant)

                progress <- function(n){
                    utils::setTxtProgressBar(pb, n)
                }

                opts <- list(progress=progress)
                start <- Sys.time()
                print(base::paste(label, " (", Sys.time(), ")", sep=""))
                pb <- txtProgressBar(max=iterations, style=3)
                res_LFQ <- foreach::foreach(i=1:nrow(LFQ_peptide_quant),
                                            .options.snow=opts) %dopar% {
                    c1 <- features$Sequence == unique_seq[i]
                    c2 <- features$Modifications == unique_mod[i]
                    sel <- which(c1 & c2)
                    testdat = features_quant[sel, ]

                    # If more than 5 peptides are available select TopN peptides
                    # according to intensity over all samples (select highest
                    # abundant peptides as quantifications are more accurate)
                    if(nrow(testdat) > TopN){
                        total_abundance <- rowSums(testdat, na.rm=TRUE)
                        loc <- order(total_abundance, decreasing=TRUE)[1:TopN]
                        testdat <- testdat[loc, ]
                    }

                    if(nrow(testdat) > 1){
                        res <- calculate_LFQ(
                            peptide_quant_data=testdat,
                            min_num_ratios=1,
                            num_ratio_samples=num_ratio_samples
                        )
                    }

                    else{
                        res <- as.numeric(testdat)
                        names(res) <- base::paste("V", 1:length(res), sep="")
                    }

                    return(res)
                }

                snow::stopCluster(cl)
                close(pb)

                # Combine results
                for(i in 1:length(res_LFQ)){
                    if(length(res_LFQ[[i]]) > 0){
                        data.table::set(LFQ_peptide_quant, as.integer(i),
                                        as.integer(1:ncol(LFQ_peptide_quant)),
                                        as.list(as.numeric(res_LFQ[[i]])))
                    }
                }

                # Add information about number of quant features per protein
                aux <- base::paste(features$Sequence[sel],
                                   features$Modifications[sel], sep="_")
                count_quant_features <- plyr::count(sort(aux))

                loc <- match(unique_peptides, count_quant_features$x)
                LFQ_peptide_quant <- base::data.frame(
                    Sequence=unique_seq,
                    Modifications=unique_mod,
                    Protein=features$Protein[match(unique_seq,
                                                   features$Sequence)],
                    num_quant_features=count_quant_features$freq[loc],
                    LFQ_peptide_quant
                )

                cnames <- base::gsub("^X", "", colnames(LFQ_peptide_quant))
                colnames(LFQ_peptide_quant) <- cnames

                end <- Sys.time()
                print(base::paste("Finished peptide LFQ-quantification (",
                                  Sys.time(), ")", sep=""))
                print(end - start)
                # Replace 0 by NA
                LFQ_peptide_quant[LFQ_peptide_quant == 0] <- NA

                return(LFQ_peptide_quant)
            }

            if(ncol(feat_w_backgr_int) <= 10){
                num_ratio_samples <- NA
            }

            if(ncol(feat_w_backgr_int) > 10){
                num_ratio_samples <- 6
            }

            LFQ_peptide_quant_with_background <- LFQ_peptide_quant_process(
                features,
                feat_w_backgr_int,
                n_cores,
                label="Perform peptide LFQ-quantification for non-imputed data",
                num_ratio_samples=num_ratio_samples
            )
            LFQ_pept_quant_w_backgrnd_imputed <- LFQ_peptide_quant_process(
                features,
                feat_w_backgr_int_imputed,
                n_cores,
                label="Perform peptide LFQ-quantification for imputed data",
                num_ratio_samples=num_ratio_samples
            )

            save(LFQ_peptide_quant_with_background,
                LFQ_pept_quant_w_backgrnd_imputed,
                file="Peptide_LFQ_temp.RData")
        }

        # Perform protein level quantification
        if(calc_protein_LFQ == TRUE){ # XXX: Having a deja vu
            # Implementation of the MaxLFQ algorithm. num_ratio_samples
            # indicates between how many samples the ratio matrices should be
            # determined. If num_ratio_samples is set to NA it will perform
            # least-square analysis between all samples. If num_ratio_samples is
            # set to a number < number of samples least square analysis is
            # performed for randomly picked n (=num_ratio_samples) samples to
            # reduced computation time
            LFQ_protein_quant_process <- function(
                features,
                features_quant,
                n_cores,
                label="Perform LFQ-quantification",
                num_ratio_samples=NA,
                TopN=5,
                seed=1
            ){
                set.seed(seed)
                # Prepare for MaxLFQ algorithm
                calculate_LFQ <- function(peptide_quant_data, min_num_ratios=2,
                                          num_ratio_samples=NA){
                    # XXX: AGAIN?
                    # Error function which is used to optimize ratios
                    least_square_error <- function(par, ratio_mat){
                        sum <- 0
                        vals <- NULL

                        if(ncol(ratio_mat) > 1){
                            for(c in 1:(ncol(ratio_mat) - 1)){
                                for(r in (c + 1):nrow(ratio_mat)){
                                    val <- (base::log2(ratio_mat[r, c])
                                            - base::log2(par[r])
                                            + base::log2(par[c]))^2
                                    vals <- append(vals, val)

                                    if(!is.na(val)){
                                        sum <- sum + val
                                    }
                                }
                            }
                        }

                        if(sum == 0){
                            sum <-  NA
                        }

                        return(sum)
                    }

                    # Unlog intensities
                    peptide_quant_data <- 2^peptide_quant_data
                    # Calculate summed intensities per sample
                    totalsum_per_sample <- colSums(peptide_quant_data,
                                                   na.rm=TRUE)
                    # Determine median ratio matrix between all samples
                    mat <- matrix(nrow=ncol(peptide_quant_data),
                                  ncol=ncol(peptide_quant_data))
                    ratio_mat <- base::as.data.frame()

                    if(ncol(ratio_mat) > 1){
                        for(c in 1:(ncol(ratio_mat) - 1)){
                            for(r in (c + 1):nrow(ratio_mat)){
                                ratios <- (peptide_quant_data[, r]
                                           / peptide_quant_data[, c])
                                loc <- which(!is.na(ratios))

                                if(length(loc) >= min_num_ratios){
                                    ratio_mat[r, c] <- stats::median(ratios,
                                                                     na.rm=TRUE)
                                }
                            }
                        }
                    }

                    # Calculate ratios over all samples
                    if(is.na(num_ratio_samples)){
                        # Define start parameter
                        start_par <- c(rep(1, ncol(ratio_mat)))
                        # Now find optimum in ratios to best recover true
                        # observed ratios between samples
                        res_ratio <- NULL
                        try(res_ratio <- stats::optim(par=start_par,
                                                      fn=least_square_error,
                                                      ratio_mat=ratio_mat,
                                                      lower=1, upper=100,
                                                      method="L-BFGS-B"),
                            silent=TRUE)

                        if(!is.null(res_ratio)){
                            # Normalize ratios to sample with highest intensity
                            mx <- max(totalsum_per_sample)
                            loc <- which(totalsum_per_sample == mx)
                            ratio_norm <- res_ratio$par / res_ratio$par[loc]
                            # Finally calculate log2 lfq protein intensities per
                            # sample
                            lfq <- base::log2(ratio_norm
                                              * totalsum_per_sample[loc])
                            # Remove quant values for samples were no ratios
                            # were available
                            c1 <- colSums(!is.na(ratio_mat)) == 0
                            c2 <- rowSums(!is.na(ratio_mat)) == 0

                            if(any(c1 & c2)){
                                sel <- as.numeric(which(c1 & c2)
                                lfq[sel] <- NA
                            }
                        }

                        else{
                            lfq <- rep(NA, ncol(peptide_quant_data))
                        }
                    }

                    # Calculate ratios only for a subset of samples
                    else{
                        mat <- matrix(nrow=ncol(peptide_quant_data),
                                      ncol=ncol(peptide_quant_data), 0)
                        lfq <- base::as.data.frame(mat)

                        for(s in 1:ncol(peptide_quant_data)){
                            # Randomly select up to 6 other samples from list of
                            # samples which also contain quantifications
                            cs <- colSums(peptide_quant_data, na.rm=TRUE)
                            aux <- which(cs > 0)
                            samples_with_quant <- as.numeric(aux)
                            loc <- which(samples_with_quant != s)
                            samples_with_quant <- samples_with_quant[loc]

                            if(length(samples_with_quant) > 0){
                                len <- length(samples_with_quant)
                                samps_comparison_x_samp <- sample(
                                    samples_with_quant,
                                    ifelse(len > num_ratio_samples,
                                           num_ratio_samples, len)
                                )
                            }

                            else{
                                samps_comparison_x_samp <- sample(
                                    c(1:ncol(peptide_quant_data))[-i],
                                    size=num_ratio_samples,
                                    replace=FALSE
                                )
                            }

                            loc <- c(s, sort(samps_comparison_x_samp))
                            ratio_mat_temp <- ratio_mat[loc, loc]
                            # To few of selected random samples show an observed
                            # intensity ratio but protein is quantified in
                            # current sample
                            c1 <- length(which(!is.na(ratio_mat_temp))) < 3
                            c2 <- any(!is.na(peptide_quant_data[, s]))

                            if(c1 & c2){
                                # Randomly select other samples
                            }

                            totalsum_per_sample_temp <- totalsum_per_sample[loc]
                            # Define start parameter
                            start_par <- c(rep(1, ncol(ratio_mat_temp)))
                            # Now find optimum in ratios to best recover true
                            # observed ratios between samples
                            res_ratio <- NULL
                            try(res_ratio <- stats::optim(
                                par=start_par,
                                fn=least_square_error,
                                ratio_mat=ratio_mat_temp,
                                lower=1,
                                upper=100,
                                method="L-BFGS-B"
                            ), silent=TRUE)

                            if(!is.null(res_ratio)){
                                # Normalize ratios to sample with highest
                                # intensity
                                mx <- max(totalsum_per_sample_temp, na.rm=TRUE)
                                loc <- which(totalsum_per_sample_temp == mx)
                                ratio_norm <- res_ratio$par / res_ratio$par[loc]
                                # Finally calculate log2 lfq protein intensities
                                # per sample
                                prod <- (ratio_norm
                                         * totalsum_per_sample_temp[loc])
                                temp_lfq <- base::log2(prod)

                                val <- c(s, sort(samps_comparison_x_samp))
                                data.table::set(lfq, as.integer(s),
                                                as.integer(val),
                                                as.list(temp_lfq))
                            }
                        }

                        lfq[lfq == 0] <- NA
                        lfq <- matrixStats::colMedians(as.matrix(lfq),
                                                       na.rm=TRUE)

                        # Remove quant values for samples were no ratios were
                        # available
                        c1 <- colSums(!is.na(ratio_mat)) == 0
                        c2 <- rowSums(!is.na(ratio_mat)) == 0

                        if(any(c1 & c2)){
                            sel <- as.numeric(which(c1 & c2))
                            lfq[sel] <- NA
                        }
                    }

                    return(lfq)
                }

                # XXX: AGAIN??
                # txtProgressBar from package pbarETA (Francesco Napolitano)
                # License: LGPL-3
                txtProgressBar <- function(min=0, max=1, initial=0, char="=",
                                           width=NA, title, label, style=3,
                                           file=""){
                    formatTime <- function(seconds){
                        if(seconds == Inf || is.nan(seconds) || is.na(seconds)){
                            return("NA")
                        }

                        seconds <- round(seconds)
                        sXmin <- 60
                        sXhr <- sXmin * 60
                        sXday <- sXhr * 24
                        sXweek <- sXday * 7
                        sXmonth <- sXweek * 4.22
                        sXyear <- sXmonth * 12
                        years <- floor(seconds / sXyear)
                        seconds <- seconds - years * sXyear
                        months <- floor(seconds / sXmonth)
                        seconds <- seconds - months * sXmonth
                        weeks <- floor(seconds / sXweek)
                        seconds <- seconds - weeks * sXweek
                        days <- floor(seconds / sXday)
                        seconds <- seconds - days * sXday
                        hours <- floor(seconds / sXhr)
                        seconds <- seconds - hours * sXhr
                        minutes <- floor(seconds / sXmin)
                        seconds <- seconds - minutes * sXmin
                        ETA <- c(years, months, days, hours, minutes, seconds)
                        startst <- which(ETA > 0)[1]

                        if(is.na(startst)){
                            startst <- 6
                        }

                        starts <- min(startst, 4)
                        fmtstr <- rep("%02d", length(ETA))[startst:length(ETA)]
                        fmtstr <- base::paste(fmtstr, collapse=":")
                        lst <- as.list(c(as.list(fmtstr),
                                         ETA[startst:length(ETA)]))

                        return(do.call(sprintf, lst))
                    }

                    c1 <- !identical(file, "")
                    c2 <- !(inherits(file, "connection") && isOpen(file))

                    if(c1 && c2){
                        stop("'file' must be \"\" or an open connection object")
                    }

                    if(!style %in% 1L:3L){
                        style <- 1
                    }

                    .val <- initial
                    .killed <- FALSE
                    .nb <- 0L
                    .pc <- -1L
                    .time0 <- NA
                    .timenow <- NA
                    .firstUpdate <- T
                    nw <- nchar(char, "w")

                    if(is.na(width)){
                        width <- getOption("width")

                        if(style == 3L){
                            width <- width - 10L
                        }

                        width <- trunc(width / nw)
                    }

                    if(max <= min){
                        stop("must have 'max' > 'min'")
                    }

                    up1 <- function(value){
                        if(!is.finite(value) || value < min || value > max){
                            return()
                        }

                        .val <<- value
                        nb <- round(width * (value - min) / (max - min))

                        if(.nb < nb){
                            cat(base::paste(rep.int(char, nb - .nb),
                                            collapse=""), file=file)
                            utils::flush.console()
                        }

                        else if(.nb > nb){
                            cat("\r", base::paste(rep.int(" ", .nb * nw),
                                                  collapse=""), "\r",
                                base::paste(rep.int(char, nb), collapse=""),
                                            sep="", file=file)
                            utils::flush.console()
                        }

                        .nb <<- nb
                    }

                    up2 <- function(value){
                        if(!is.finite(value) || value < min || value > max){
                            return()
                        }

                        .val <<- value
                        nb <- round(width * (value - min) / (max - min))

                        if(.nb <= nb){
                            cat("\r", base::paste(rep.int(char, nb),
                                                  collapse=""), sep="",
                                file=file)
                            utils::flush.console()
                        }

                        else{
                            cat("\r", base::paste(rep.int(" ", .nb * nw),
                                                  collapse=""), "\r",
                                base::paste(rep.int(char, nb), collapse=""),
                                            sep="", file=file)
                            utils::flush.console()
                        }

                        .nb <<- nb
                    }

                    up3 <- function(value, calledOnCreation=FALSE){
                        timenow <- proc.time()[["elapsed"]]

                        if(!calledOnCreation && .firstUpdate){
                            .time0 <<- timenow
                            .timenow <<- timenow
                            .firstUpdate <<- FALSE
                        }

                        if(!is.finite(value) || value < min || value > max){
                            return()
                        }

                        .val <<- value
                        nb <- round(width * (value - min) / (max - min))
                        pc <- round(100 * (value - min) / (max - min))

                        if (nb == .nb && pc == .pc && timenow - .timenow < 1){
                            return()
                        }

                        .timenow <<- timenow
                        span <- timenow - .time0
                        timeXiter <- span / (.val - min)
                        ETA <- (max - .val) * timeXiter
                        ETAstr <- formatTime(ETA)
                        cat(base::paste(c("\r  |",
                                          rep.int(" ", nw * width + 6)),
                                        collapse = ""), file=file)
                        cat(base::paste(c("\r  |", rep.int(char, nb),
                                          rep.int(" ", nw * (width - nb)),
                                          sprintf("| %3d%%", pc), ", ETA ",
                                          ETAstr), collapse=""), file=file)
                        utils::flush.console()
                        .nb <<- nb
                        .pc <<- pc
                    }

                    getVal <- function(){
                        .val
                    }

                    kill <- function(){
                        if (!.killed) {
                            cat("\n", file=file)
                            utils::flush.console()
                            .killed <<- TRUE
                        }
                    }

                    up <- switch(style, up1, up2, up3)
                    up(initial, TRUE)
                    structure(list(getVal=getVal, up=up, kill=kill),
                              class="txtProgressBar")
                }

                loc <- which(!grepl(";|,", features$Protein))
                unique_proteins <- sort(unique(features$Protein[loc]))
                unique_proteins <- unique_proteins[which(unique_proteins != "")]

                mat <- matrix(ncol=ncol(features_quant),
                              nrow=length(unique_proteins), 0)
                LFQ_protein_quant <- base::as.data.frame(mat)
                colnames(LFQ_protein_quant) <- colnames(features_quant)
                rownames(LFQ_protein_quant) <- unique_proteins

                # Perform LFQ quantification
                cl <- snow::makeCluster(n_cores)
                doSNOW::registerDoSNOW(cl)
                iterations <- nrow(LFQ_protein_quant)

                progress <- function(n){
                    utils::setTxtProgressBar(pb, n)
                }

                opts <- list(progress=progress)
                start <- Sys.time()
                print(base::paste(label, " (", Sys.time(), ")", sep=""))
                pb <- txtProgressBar(max=iterations, style=3)
                res_LFQ <- foreach::foreach(i=1:nrow(LFQ_protein_quant),
                                            .options.snow=opts) %dopar% {
                    sel <- which(features$Protein == unique_proteins[i])
                    testdat = features_quant[sel, ]

                    # If more than 5 peptides are available select TopN peptides
                    # according to intensity over all samples (select highest
                    # abundant peptides as quantifications are more accurate)
                    if(nrow(testdat) > TopN){
                        total_abundance <- rowSums(testdat, na.rm=TRUE)
                        testdat <- testdat[order(total_abundance,
                                                 decreasing=TRUE)[1:TopN], ]
                    }

                    if(nrow(testdat) >= 2){
                        res <- calculate_LFQ(
                            peptide_quant_data=testdat,
                            min_num_ratios=2,
                            num_ratio_samples=num_ratio_samples
                        )
                    }

                    else{
                        res <- as.numeric(rep(NA, ncol(features_quant)))
                        names(res) <- base::paste("V", 1:length(res), sep="")
                    }

                    return(res)
                }

                snow::stopCluster(cl)
                close(pb)

                #Combine results
                for(i in 1:length(res_LFQ)){
                    if(length(res_LFQ[[i]]) > 0){
                        data.table::set(LFQ_protein_quant, as.integer(i),
                                        as.integer(1:ncol(LFQ_protein_quant)),
                                        s.list(as.numeric(res_LFQ[[i]])))
                    }
                }

                # Add information about number of quant features per protein
                loc <- which(!grepl(";|\\||,", features$Protein))
                count_quant_features <- plyr::count(features$Protein[loc])

                loc <- match(rownames(LFQ_protein_quant),
                             count_quant_features$x)
                LFQ_protein_quant <- base::data.frame(
                    num_quant_features=count_quant_features$freq[loc],
                    LFQ_protein_quant
                )

                end <- Sys.time()
                print(base::paste("Finished LFQ-quantification (", Sys.time(),
                                  ")", sep=""))
                print(end - start)
                # Replace 0 by NA
                LFQ_protein_quant[LFQ_protein_quant == 0] <- NA

                return(LFQ_protein_quant)
            }

            if(ncol(feat_w_backgr_int) <= 10){
                num_ratio_samples <- NA
            }

            if(ncol(feat_w_backgr_int) > 10){
                num_ratio_samples <- 6
            }

            if(calc_peptide_LFQ == FALSE){
                LFQ_quant_with_background <- LFQ_protein_quant_process(
                    features,
                    feat_w_backgr_int,
                    n_cores,
                    label="Perform LFQ-quantification for non-imputed data",
                    num_ratio_samples=num_ratio_samples
                )
                LFQ_quant_with_background_imputed <- LFQ_protein_quant_process(
                    features,
                    feat_w_backgr_int_imputed,
                    n_cores,
                    label="Perform LFQ-quantification for imputed data",
                    num_ratio_samples=num_ratio_samples
                )
            }

            # Use peptide LFQ for calcualting protein LFQ
            else{
                loc1 <- c(5:ncol(LFQ_peptide_quant_with_background))
                loc2 <- c(5:ncol(LFQ_pept_quant_w_backgrnd_imputed))
                LFQ_quant_with_background <- LFQ_protein_quant_process(
                    LFQ_peptide_quant_with_background[, c(1:4)],
                    LFQ_peptide_quant_with_background[, loc1],
                    n_cores,
                    label="Perform LFQ-quantification for non-imputed data",
                    num_ratio_samples=num_ratio_samples
                )
                LFQ_quant_with_background_imputed <- LFQ_protein_quant_process(
                    LFQ_pept_quant_w_backgrnd_imputed[, c(1:4)],
                    LFQ_pept_quant_w_backgrnd_imputed[, loc2],
                    n_cores,
                    label="Perform LFQ-quantification for imputed data",
                    num_ratio_samples=num_ratio_samples
                )
            }
        }

        # Perform Top3 and Total quantification
        cl <- parallel::makeCluster(ifelse(n_cores < 4, n_cores, 4))
        doParallel::registerDoParallel(cl)

        res <- foreach::foreach(i=1:4) %dopar% {
            if(i == 1){
                # Perform Top3 protein quantification
                res <- Top3_Protein_Quant(
                    features=features,
                    feat_sample_mat_requant=feat_w_backgr_int,
                    Quant_pvals=pval_sig_w_bckgrnd_quant,
                    S2B=S2B,
                    Alignment_scores=align_scores_peaks_correct,
                    use_overlapping=TRUE,
                    min_peps=1
                )
            }

            if(i == 2){
                # Perform Total protein quantification
                res <- Total_Protein_Quant(
                    features=features,
                    feat_sample_mat_requant=feat_w_backgr_int,
                    Quant_pvals=pval_sig_w_bckgrnd_quant,
                    S2B=S2B,
                    Alignment_scores=align_scores_peaks_correct,
                    use_overlapping=TRUE,
                    min_peps=1
                )
            }

            if(i == 3){
                # Perform Top3 protein quantification for imputed data
                res <- Top3_Protein_Quant(
                    features=features,
                    feat_sample_mat_requant=feat_w_backgr_int_imputed,
                    Quant_pvals=pval_sig_w_bckgrnd_quant,
                    Alignment_scores=align_scores_peaks_correct,
                    S2B=S2B,
                    use_overlapping=TRUE,
                    min_peps=1
                )
            }

            if(i == 4){
                # Perform Total protein quantification for imputed data
                res <- Total_Protein_Quant(
                    features=features,
                    feat_sample_mat_requant=feat_w_backgr_int_imputed,
                    Quant_pvals=pval_sig_w_bckgrnd_quant,
                    Alignment_scores=align_scores_peaks_correct,
                    S2B=S2B,
                    use_overlapping=T,
                    min_peps=1
                )
            }

            return(res)
        }

        parallel::stopCluster(cl)

        Top3_quant_with_background <- res[[1]]
        Total_quant_with_background <- res[[2]]
        Top3_quant_with_background_imputed <- res[[3]]
        Total_quant_with_background_imputed <- res[[4]]

        # Match Gene names to Uniprot Identifier
        c1 <- any(colnames(MaxQ_protein_groups) == "Gene.names")
        c2 <- any(colnames(MaxQ_protein_groups) == "Protein.IDs")

        if(c1 & c2){
            temp <- MaxQ_protein_groups[, c("Gene.names", "Protein.IDs")]
            temp <- temp[which(temp$Gene.names != ""), ]

            upid <- unique(as.character(stringr::str_split(temp$Protein.IDs,
                                                           ";", simplify=TRUE)))
            UniProt_to_GeneName <- base::data.frame(UniProt_ID=upid,
                                                    Gene_Name="")
            gname <- as.character(UniProt_to_GeneName$Gene_Name)
            UniProt_to_GeneName$Gene_Name <- gname

            for(i in 1:nrow(UniProt_to_GeneName)){
                aux <- grepl(UniProt_to_GeneName$UniProt_ID[i],
                             temp$Protein.IDs)
                gn <- temp$Gene.names[which(aux)]
                UniProt_to_GeneName$Gene_Name[i] <- as.character(gn)
            }

            loc <- match(rownames(Top3_quant_with_background),
                         UniProt_to_GeneName$UniProt_ID)
            Top3_quant_with_background <- base::data.frame(
                Gene_Name=UniProt_to_GeneName$Gene_Name[loc],
                UniProt_Identifier=rownames(Top3_quant_with_background),
                Top3_quant_with_background
            )
            rownames(Top3_quant_with_background) <- c()

            loc <- match(rownames(Total_quant_with_background),
                         UniProt_to_GeneName$UniProt_ID)
            Total_quant_with_background <- base::data.frame(
                Gene_Name=UniProt_to_GeneName$Gene_Name[loc],
                UniProt_Identifier=rownames(Total_quant_with_background),
                Total_quant_with_background
            )
            rownames(Total_quant_with_background) <- c()

            loc <- match(rownames(Top3_quant_with_background_imputed),
                         UniProt_to_GeneName$UniProt_ID)
            Top3_quant_with_background_imputed <- base::data.frame(
                Gene_Name=UniProt_to_GeneName$Gene_Name[loc],
                UniProt_Identifier=rownames(Top3_quant_with_background_imputed),
                Top3_quant_with_background_imputed
            )
            rownames(Top3_quant_with_background_imputed) <- c()
            rn <- rownames(Total_quant_with_background_imputed)
            loc <- match(rn, UniProt_to_GeneName$UniProt_ID)

            Total_quant_with_background_imputed <- base::data.frame(
                Gene_Name=UniProt_to_GeneName$Gene_Name[loc],
                UniProt_Identifier=rn,
                Total_quant_with_background_imputed
            )
            rownames(Total_quant_with_background_imputed) <- c()

            if(calc_peptide_LFQ == TRUE){
                upid <- LFQ_peptide_quant_with_background$Protein
                loc1 <- match(upid, UniProt_to_GeneName$UniProt_ID)
                loc2 <- c(1, 2, 4, 5:ncol(LFQ_peptide_quant_with_background))
                LFQ_peptide_quant_with_background <- base::data.frame(
                    Gene_Name=UniProt_to_GeneName$Gene_Name[loc1],
                    UniProt_Identifier=upid,
                    LFQ_peptide_quant_with_background[, loc2]
                )
                rownames(LFQ_peptide_quant_with_background) <- c()
                colnames(LFQ_peptide_quant_with_background) <- base::gsub(
                    "^X",
                    "",
                    colnames(LFQ_peptide_quant_with_background)
                )

                loc1 <- match(LFQ_pept_quant_w_backgrnd_imputed$Protein,
                              UniProt_to_GeneName$UniProt_ID)
                loc2 <- c(1, 2, 4,
                          5:ncol(LFQ_pept_quant_w_backgrnd_imputed))
                aux <- LFQ_pept_quant_w_backgrnd_imputed$Protein
                LFQ_pept_quant_w_backgrnd_imputed <- base::data.frame(
                    Gene_Name=UniProt_to_GeneName$Gene_Name[loc1],
                    UniProt_Identifier=aux,
                    LFQ_pept_quant_w_backgrnd_imputed[, loc2]
                )
                rownames(LFQ_pept_quant_w_backgrnd_imputed) <- c()
                colnames(LFQ_pept_quant_w_backgrnd_imputed) <- base::gsub(
                    "^X",
                    "",
                    colnames(LFQ_pept_quant_w_backgrnd_imputed)
                )
            }

            if(calc_protein_LFQ == TRUE){
                loc <- match(rownames(LFQ_quant_with_background),
                             UniProt_to_GeneName$UniProt_ID)
                LFQ_quant_with_background <- base::data.frame(
                    Gene_Name=UniProt_to_GeneName$Gene_Name[loc],
                    UniProt_Identifier=rownames(LFQ_quant_with_background),
                    LFQ_quant_with_background
                )
                rownames(LFQ_quant_with_background) <- c()
                colnames(LFQ_quant_with_background) <- base::gsub(
                    "^X",
                    "",
                    colnames(LFQ_quant_with_background)
                )

                upid <- rownames(LFQ_quant_with_background_imputed)
                loc <- match(upid, UniProt_to_GeneName$UniProt_ID)
                LFQ_quant_with_background_imputed <- base::data.frame(
                    Gene_Name=UniProt_to_GeneName$Gene_Name[loc],
                    UniProt_Identifier=upid,
                    LFQ_quant_with_background_imputed
                )
                rownames(LFQ_quant_with_background_imputed) <- c()
                colnames(LFQ_quant_with_background_imputed) <- base::gsub(
                    "^X",
                    "",
                    colnames(LFQ_quant_with_background_imputed)
                )
            }

        }

        if(calc_protein_LFQ == TRUE){
            fname <- base::paste(path_to_features,
                                 "/Proteins_quantification_LFQ",
                                 output_file_names_add, ".tab", sep="")
            utils::write.table(x=LFQ_quant_with_background, file=fname,
                               row.names=FALSE, sep="\t")
            fname <- base::paste(path_to_features,
                                 "/Proteins_quantification_LFQ_imputed",
                                 output_file_names_add, ".tab", sep="")
            utils::write.table(x=LFQ_quant_with_background_imputed, file=fname,
                               row.names=FALSE, sep="\t")
        }

        fname <- base::paste(path_to_features, "/Proteins_quantification_Top3",
                             output_file_names_add, ".tab", sep="")
        utils::write.table(x=Top3_quant_with_background, file=fname,
                           row.names=FALSE, sep="\t")
        fname <- base::paste(path_to_features,
                             "/Proteins_quantification_Top3_imputed",
                             output_file_names_add, ".tab", sep="")
        utils::write.table(x=Top3_quant_with_background_imputed, file=fname,
                           row.names=FALSE, sep="\t")
        fname <- base::paste(path_to_features, "/Proteins_quantification_Total",
                             output_file_names_add, ".tab", sep="")
        utils::write.table(x=Total_quant_with_background, file=fname,
                           row.names=FALSE, sep="\t")
        fname <- base::paste(path_to_features,
                             "/Proteins_quantification_Total_imputed",
                             output_file_names_add, ".tab", sep="")
        utils::write.table(x=Total_quant_with_background_imputed, file=fname,
                           row.names=FALSE, sep="\t")

        if(calc_peptide_LFQ == TRUE){
            fname <- base::paste(path_to_features,
                                 "/Peptides_quantification_LFQ",
                                 output_file_names_add, ".tab", sep="")
            utils::write.table(x=LFQ_peptide_quant_with_background, file=fname,
                               row.names=FALSE, sep="\t")
            fname <- base::paste(path_to_features,
                                 "/Peptides_quantification_LFQ_imputed",
                                 output_file_names_add, ".tab", sep="")
            utils::write.table(x=LFQ_pept_quant_w_backgrnd_imputed,
                               file=fname, row.names=FALSE, sep="\t")
        }

        print(paste0(Sys.time(), " Protein-level aggregation finished"))
    }

    setwd(path_to_features)

    # Finally save feature level quantification
    print(paste0(Sys.time(), " Store all results"))
    fname <- base::paste(path_to_features, "/Features", output_file_names_add,
                         ".tab", sep="")
    utils::write.table(x=features, file=fname, row.names=FALSE, sep="\t")
    fname <- base::paste(path_to_features, "/Features_quantification",
                         output_file_names_add, ".tab", sep="")
    utils::write.table(x=feat_w_backgr_int, file=fname, row.names=TRUE,
                       sep="\t")
    fname <- base::paste(path_to_features, "/Features_quantification_imputed",
                         output_file_names_add, ".tab", sep="")
    utils::write.table(x=feat_w_backgr_int_imputed, file=fname, row.names=TRUE,
                       sep="\t")
    fname <- base::paste(path_to_features, "/Features_quantification_pvals",
                         output_file_names_add, ".tab", sep="")
    utils::write.table(x=pval_sig_w_bckgrnd_quant, file=fname, row.names=TRUE,
                       sep="\t")
    fname <- base::paste(path_to_features, "/Features_quantification_ioncount",
                         output_file_names_add, ".tab", sep="")
    utils::write.table(x=Icount_feat_w_bkgr_int, file=fname, row.names=TRUE,
                       sep="\t")
    fname <- base::paste(path_to_features, "/Features_quantification_S2B",
                         output_file_names_add, ".tab", sep="")
    utils::write.table(x=S2B, file=fname, row.names=TRUE, sep="\t")
    fname <- base::paste(path_to_features,
                         "/Features_quantification_variability_score",
                         output_file_names_add, ".tab", sep="")
    utils::write.table(x=align_var_score, file=fname, row.names=TRUE, sep="\t")
    fname <- base::paste(path_to_features,
                         "/Features_quantification_alignment_score",
                         output_file_names_add, ".tab", sep="")
    utils::write.table(x=align_scores_peaks_correct, file=fname, row.names=TRUE,
                       sep="\t")
    fname <- base::paste(path_to_features,
                         "/Features_quantification_mono_iso_alignment_score",
                         output_file_names_add, ".tab", sep="")
    utils::write.table(x=mono_iso_alignment_summary, file=fname, row.names=TRUE,
                       sep="\t")

    save(QC_data, file="Temporary_files/Feature_quantification_QC_data.RData")
    options(warn=0)
}

#' Peak decision algorithm only designed for internal use
#' @param features_select Features on which peak-selection should be performed
#' @param peak_quant List of determined peaks per feature
#' @param samples String vector of all samples
#' @param s Integer indicating index of current sample
#' @param RT_correction_factors RT corrections
#' @param mz_correction_factors mz corrections
#' @param features_intensity_sample features_intensity_sample
#' @param Ioncount_sample Ioncount_sample
#' @param feat_w_bkgrnd_inty_samp
#' feat_w_bkgrnd_inty_samp
#' @param Icount_w_bkgr_smp Icount_w_bkgr_smp
#' @param peak_selected_sample peak_selected_sample
#' @param delta_mz delta_mz
#' @param delta_rt delta_rt
#' @param peak_min_ion_count peak_min_ion_count
#' @param chunk chunk
#' @param num_chunks num_chunks
#' @param progress Progressbar
#' @import data.table
#' @export
#' @details Peak decision algorithm
peak_decision <- function(features_select, peak_quant, samples, s,
                          RT_correction_factors, mz_correction_factors,
                          features_intensity_sample, Ioncount_sample,
                          feat_w_bkgrnd_inty_samp,
                          Icount_w_bkgr_smp, peak_selected_sample,
                          delta_mz, delta_rt, peak_min_ion_count, chunk=NULL,
                          num_chunks=NULL, progress=TRUE){
    # Prevent issues during R CMD check
    ..s <- s
    rm(..s)
    ..known_peaks_indices <- 1
    rm(..known_peaks_indices)
    ..o <- 1
    rm(..o)

    if(progress == TRUE){
        title <- base::paste("Prepare for peak selection -", samples[s])
        label <- base::paste(round(0 / 1 * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title=title, label=label, min=0, max=1,
                                   width=300)
        close(pb)
    }

    # Next, decide for all other feature quantifications for which at least 1
    # peak was available which peak quantification should be used
    max <- nrow(features_select)

    if(progress == TRUE){
        if(!is.null(chunk) & !is.null(num_chunks)){
            title <- base::paste("Select peaks (", chunk, "/", num_chunks,
                                 ") - ", samples[s], sep="")
            label <- base::paste(round(0 / max * 100, 0), "% done")
            pb <- tcltk::tkProgressBar(title=title, label=label, min=0, max=max,
                                       width=300)
        }

        else{
            title <- base::paste("Select peaks -", samples[s])
            label <- base::paste(round(0 / max * 100, 0), "% done")
            pb <- tcltk::tkProgressBar(title=title, label=label, min=0, max=max,
                                       width=300)
        }
    }

    start_time <- Sys.time()
    updatecounter <- 0
    time_require <- 0

    for(i in 1:nrow(features_select)){
        c1 <- !is.na(peak_quant$Standard$num_peaks_with_background[i, ..s])
        c2 <- as.numeric(peak_quant$Standard$num_peaks_with_background[i, ..s])
        c3_1 <- !is.na(peak_quant$Peak_1$num_peaks_with_background[i, ..s])
        c3_2 <- as.numeric(peak_quant$Peak_1$num_peaks_with_background[i, ..s])
        ndp <- ifelse(c1, c2, ifelse(c3_1, c3_2, NA))

        if(!is.na(ndp)){ # Peaks detected?
            # Check if more than 1 peak was detected and true peak is unknown
            if(ndp > 1){
                # More than 1 peaks but true peak not known
                if(peak_quant$Peak_1$correct_peak_w_background[i,..s] == 0){
                    # Check if for any sample the correct peak is known
                    c1 <- peak_quant$Peak_1$correct_peak_w_background[i, ]
                    num_known_peaks <- length(which(c1 != 0))

                    if(num_known_peaks > 0){
                        # Now check where correct known peak usually is located
                        # in samples were true location is known
                        pk <- peak_quant$Peak_1
                        c1 <- pk$correct_peak_w_background[i, ] == 1
                        kpi <- which(c1)

                        RT_in_known <- pk$Peak_rt_with_background[i, ..kpi]
                        RT_in_known <- as.numeric(RT_in_known)
                        mz_in_known <- pk$Peak_mz_with_background[i, ..kpi]
                        mz_in_known <- as.numeric(mz_in_known)

                        # Get RT and mz for detected peaks in current sample
                        detected_peaks_rt <- vector("numeric", as.numeric(ndp))
                        detected_peaks_mz <- vector("numeric", as.numeric(ndp))

                        names(detected_peaks_rt) <- 1:as.numeric(ndp)
                        names(detected_peaks_mz) <- 1:as.numeric(ndp)

                        for(p in 1:length(detected_peaks_rt)){
                            pk <- peak_quant[[p]]
                            dprt <- pk$Peak_rt_with_background[i, ..s]
                            detected_peaks_rt[p] <- as.numeric(dprt)
                            dpmz <- pk$Peak_mz_with_background[i, ..s]
                            detected_peaks_mz[p] <- as.numeric(dpmz)
                        }

                        # NOTE: Nesting levels are getting ridiculous
                        # TODO: Ignoring line width limit until proper refactor.
                        if(length(which(!is.na(detected_peaks_rt)
                                        & !is.na(detected_peaks_mz))) > 0){
                            detected_peaks_rt <- detected_peaks_rt[which(!is.na(detected_peaks_rt))]
                            detected_peaks_mz <- detected_peaks_mz[which(!is.na(detected_peaks_mz))]

                            # Now compare observed deviations in RT and mz for
                            # known peaks with the deviation for the potential
                            # peaks in current sample. For this comparison we
                            # have to correct detected peaks RT and m/z with the
                            # RT and mz correction factors per corresponding
                            # sample. Otherwise it would be possible that
                            # current sample shows e.g. a global RT shift which
                            # could result in a major deviation from all other
                            # known samples. Thus a peak wouldn't be selected
                            # although it e.g. perfectly lies close to the
                            # expected window

                            RT_in_known_corrected <- as.numeric(RT_in_known - RT_correction_factors[i, kpi])
                            mz_in_known_corrected <- as.numeric(mz_in_known - mz_correction_factors[i, kpi])

                            detected_peaks_rt_corrected <- as.numeric(detected_peaks_rt - RT_correction_factors[i, s])
                            detected_peaks_mz_corrected <- as.numeric(detected_peaks_mz - mz_correction_factors[i, s])
                            names(detected_peaks_rt_corrected) <- names(detected_peaks_rt)
                            names(detected_peaks_mz_corrected) <- names(detected_peaks_mz)
                            # Check if the detected RT and mz of peaks are
                            # likely to be belonging to the same population
                            # (compared to known peaks). 99 % of observed data
                            # points lie within 2*sd range --> if delta RT or
                            # delta mZ > 2*sd then most likely this is the wrong
                            # peak. So we assume all peaks within this RT and mz
                            # deviation to be valid peak candidates. If these
                            # deviations are smaller then the expected RT and mz
                            # window, we expand these acceptance criteria
                            # accordingly
                            delta_RT_cut <- 3 * stats::sd(RT_in_known_corrected, na.rm=TRUE)

                            if(is.na(delta_RT_cut) | delta_RT_cut < 2 * delta_rt){
                                delta_RT_cut <- 2*delta_rt
                            }

                            delta_mz_cut <- 3 * stats::sd(mz_in_known_corrected, na.rm=TRUE)

                            if(is.na(delta_mz_cut) | delta_mz_cut < 3 * delta_mz){
                                delta_mz_cut <- 3 * delta_mz
                            }

                            dif1 <- abs(detected_peaks_rt_corrected - mean(RT_in_known_corrected, na.rm=TRUE))
                            dif2 <- abs(detected_peaks_mz_corrected - mean(mz_in_known_corrected, na.rm=TRUE))
                            within_range <- ifelse(dif1 <= delta_RT_cut & dif2 <= delta_mz_cut, TRUE, FALSE)

                            if(any(within_range)){
                                detected_peaks_rt <- detected_peaks_rt[within_range]
                                detected_peaks_mz <- detected_peaks_mz[within_range]
                                detected_peaks_rt_corrected <- detected_peaks_rt_corrected[within_range]
                                detected_peaks_mz_corrected <- detected_peaks_mz_corrected[within_range]

                                # Now check that at the position of the
                                # potential peak no other peak is present in
                                # samples with known peak
                                count_other_peaks <- sum(peak_quant$Peak_1$num_peaks_with_background[i, ..kpi] - 1)

                                if(count_other_peaks > 0){
                                    other_peaks <- base::as.data.frame(matrix(ncol=6, nrow=count_other_peaks, 0))
                                    counter <- 0

                                    # Collect peak data for all additional peaks
                                    # in samples were peak was known
                                    for(o in kpi){
                                        cur_peak_count <- as.numeric(peak_quant$Peak_1$num_peaks_with_background[i, ..o] - 1)

                                        if(cur_peak_count > 0){
                                            for(p in 1:cur_peak_count){
                                                counter <- counter + 1
                                                temp <- c(as.numeric(peak_quant[[p + 1]]$Peak_mz_with_background[i, ..o] - mz_correction_factors[i, o]),
                                                          as.numeric(peak_quant[[p + 1]]$Peak_rt_with_background[i, ..o] - RT_correction_factors[i, o]),
                                                          as.numeric(peak_quant[[p + 1]]$feat_w_backgr_int[i, ..o]),
                                                          as.numeric(peak_quant[[p + 1]]$Icount_feat_w_bkgr_int[i, ..o]))
                                                # Add distance in RT and mz to
                                                # known peak
                                                temp[5:6] <- c(abs(temp[1] - as.numeric(peak_quant[[1]]$Peak_mz_with_background[i, ..o] - mz_correction_factors[i, o])),
                                                               abs(temp[2] - as.numeric(peak_quant[[1]]$Peak_rt_with_background[i, ..o] - RT_correction_factors[i, o])))

                                                data.table::set(other_peaks, as.integer(counter), as.integer(1:6), value=list(temp[1], temp[2], temp[3], temp[4], temp[5], temp[6]))
                                            }
                                        }
                                    }

                                    # Disregard other peaks below significance
                                    # threshold
                                    other_peaks <- other_peaks[which(other_peaks$V4 > peak_min_ion_count), ]
                                    # Disregard other peaks which are too close
                                    # to the expected peak
                                    # Cutoff RT: > 1.1 * half peak width -> half RT extraction window
                                    # Cutoff mz: > 1.1 * delta_mz -> half mz extraction window
                                    other_peaks <- other_peaks[which(other_peaks$V5 > delta_mz * 1.1 | other_peaks$V6 > ((features_select$RT_length[i] / 2) * 1.1)), ]
                                    # Now check that selected peaks are not
                                    # close to other peaks in sampels with known
                                    # peak
                                    overlap <- vector("logical", length(detected_peaks_rt))

                                    for(o in 1:length(detected_peaks_rt)){
                                        overlap[o] <- ifelse(any(abs(detected_peaks_rt_corrected[o] - other_peaks$V2) <= (delta_rt / 2) & abs(detected_peaks_mz_corrected[o] - other_peaks$V1) <= delta_mz / 2, na.rm=TRUE), TRUE, FALSE)
                                    }
                                }

                                else{
                                    overlap <- vector("logical", length(detected_peaks_rt))
                                }

                                # Not overlapping with any other peak
                                if(any(overlap == FALSE)){
                                    detected_peaks_rt <- detected_peaks_rt_corrected[which(overlap == FALSE)]
                                    detected_peaks_mz <- detected_peaks_mz_corrected[which(overlap == FALSE)]
                                    RT_in_known_other_peaks <- as.numeric(peak_quant$Peak_1$Peak_rt_with_background[i, ..kpi])
                                    mz_in_known_other_peaks <- as.numeric(peak_quant$Peak_1$Peak_mz_with_background[i, ..kpi])

                                    # Select the peak which is closest to all other known peaks
                                    sum_delta_to_known_peaks_rt <- vector("numeric", length(detected_peaks_rt))
                                    sum_delta_to_known_peaks_mz <- vector("numeric", length(detected_peaks_rt))
                                    sum_delta_to_known_peaks <- vector("numeric", length(detected_peaks_rt))

                                    for(p in 1:length(detected_peaks_rt)){
                                        sum_delta_to_known_peaks_rt[p] <- sum(detected_peaks_rt[p] - RT_in_known_corrected, na.rm=TRUE) / length(RT_in_known_corrected)
                                        sum_delta_to_known_peaks_mz[p] <- sum((detected_peaks_mz[p] - mz_in_known_corrected) * 500, na.rm=TRUE) / length(RT_in_known_corrected)
                                        sum_delta_to_known_peaks[p] <- abs(sum_delta_to_known_peaks_rt[p]) + abs(sum_delta_to_known_peaks_mz[p])
                                    }

                                    selected_peak <- which(!is.na(sum_delta_to_known_peaks) & sum_delta_to_known_peaks == min(sum_delta_to_known_peaks[which(!is.na(sum_delta_to_known_peaks))], na.rm=TRUE))[1]
                                    selected_peak <- as.numeric(names(detected_peaks_rt)[selected_peak])
                                    data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant[[selected_peak]]$features_intensity[i, ..s]))
                                    data.table::set(Ioncount_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant[[selected_peak]]$Icount_feat_sample_mat[i, ..s]))
                                    data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant[[selected_peak]]$feat_w_backgr_int[i, ..s]))
                                    data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant[[selected_peak]]$Icount_feat_w_bkgr_int[i, ..s]))
                                    data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(selected_peak))
                                }

                                # None of the peaks is not overlapping with an
                                # other peak detected in other samples where
                                # correct peak was known
                                else{
                                    data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$features_intensity[i, ..s]))
                                    data.table::set(Ioncount_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_sample_mat[i, ..s]))
                                    data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$feat_w_backgr_int[i, ..s]))
                                    data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_w_bkgr_int[i, ..s]))
                                    data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(length(peak_quant)))
                                }
                            }

                            # None of the peaks is within the accepted range
                            # thus use standard window
                            else{
                                data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$features_intensity[i, ..s]))
                                data.table::set(Ioncount_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_sample_mat[i, ..s]))
                                data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$feat_w_backgr_int[i, ..s]))
                                data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_w_bkgr_int[i, ..s]))
                                data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(length(peak_quant)))
                            }
                        }

                        # None of the peaks is valid
                        else{
                            data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$features_intensity[i, ..s]))
                            data.table::set(Ioncount_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_sample_mat[i, ..s]))
                            data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$feat_w_backgr_int[i, ..s]))
                            data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_w_bkgr_int[i, ..s]))
                            data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(length(peak_quant)))
                        }
                    }

                    # In no other sample the correct peak is known thus take
                    # total window
                    else{
                        pk <- peak_quant$Standard
                        v <- pk$features_intensity[i, ..s]
                        data.table::set(features_intensity_sample,
                                        i=as.integer(i),  j=as.integer(1),
                                        value=as.numeric(v))
                        v <- pk$Icount_feat_sample_mat[i, ..s]
                        data.table::set(Ioncount_sample,
                                        i=as.integer(i), j=as.integer(1),
                                        value=as.numeric(v))
                        v <- pk$feat_w_backgr_int[i, ..s]
                        data.table::set(feat_w_bkgrnd_inty_samp,
                                        i=as.integer(i), j=as.integer(1),
                                        value=as.numeric(v))
                        v <- pk$Icount_feat_w_bkgr_int[i, ..s]
                        data.table::set(Icount_w_bkgr_smp,
                                        i=as.integer(i), j=as.integer(1),
                                        value=as.numeric(v))
                        data.table::set(peak_selected_sample,
                                        i=as.integer(i), j=as.integer(1),
                                        value=as.numeric(length(peak_quant)))
                    }
                }
            }

            # One peak was detected
            else if(ndp == 1){
                # Check if correct peak is unknown
                if(peak_quant$Peak_1$correct_peak_w_background[i, ..s] == 0){
                    # If this is the case check if detected peak is comparable
                    # to correct peaks in other samples where correct peak was
                    # known
                    # Check if for any sample the correct peak is known
                    pk <- peak_quant$Peak_1
                    cpkwb <- pk$correct_peak_w_background[i,]
                    num_known_peaks <- length(which(cpkwb != 0))

                    if(num_known_peaks > 0){
                        # Now check where correct known peak usually is located
                        # in samples were true location is known
                        kpi <- which(cpkwb == 1)

                        pkrt <- pk$Peak_rt_with_background[i, ..kpi]
                        pkmz <- pk$Peak_mz_with_background[i, ..kpi]
                        RT_in_known <- as.numeric(pkrt)
                        mz_in_known <- as.numeric(pkmz)

                        pkq_rt <- peak_quant[[1]]$Peak_rt_with_background
                        pkq_mz <- peak_quant[[1]]$Peak_mz_with_background
                        detected_peaks_rt <- as.numeric(pkq_rt[i, ..s])
                        detected_peaks_mz <- as.numeric(pkq_mz[i, ..s])

                        # NOTE: Nesting levels are getting ridiculous
                        # TODO: Ignoring line width limit until proper refactor.
                        if(length(which(!is.na(detected_peaks_rt) & !is.na(detected_peaks_mz))) > 0){
                            detected_peaks_rt <- detected_peaks_rt[which(!is.na(detected_peaks_rt))]
                            detected_peaks_mz <- detected_peaks_mz[which(!is.na(detected_peaks_mz))]

                            # Now compare observed deviations in RT and mz for
                            # known peaks with the deviation for the potential
                            # peaks in current sample. For this comparison we
                            # have to correct detected peaks RT and m/z with the
                            # RT and mz correction factors per corresponding
                            # sample. Otherwise it would be possible that
                            # current sample shows e.g. a global RT shift which
                            # could result in a major deviation from all other
                            # known samples. Thus a peak wouldn't be selected
                            # although it e.g. perfectly lies close to the
                            # expected window
                            RT_in_known_corrected <- as.numeric(RT_in_known - RT_correction_factors[i, kpi])
                            mz_in_known_corrected <- as.numeric(mz_in_known - mz_correction_factors[i, kpi])

                            detected_peaks_rt_corrected <- as.numeric(detected_peaks_rt - RT_correction_factors[i, s])
                            detected_peaks_mz_corrected <- as.numeric(detected_peaks_mz - mz_correction_factors[i, s])

                            # Check if the detected RT and mz of peaks are
                            # likely to be belonging to the same population
                            # (compared to known peaks). 99 % of observed data
                            # points lie within 3*sd range --> if delta RT or
                            # delta mZ > 2*sd then most likely this is the wrong
                            # peak. So we assume all peaks within this RT and mz
                            # deviation to be valid peak candidates. If these
                            # deviations are smaller then the expected RT and mz
                            # window, we expand these acceptance criteria
                            # accordingly
                            delta_RT_cut <- 3 * stats::sd(RT_in_known, na.rm=TRUE)

                            if(is.na(delta_RT_cut) | delta_RT_cut < 2 * delta_rt){
                                delta_RT_cut <- 2 * delta_rt
                            }

                            delta_mz_cut <- 3 * stats::sd(mz_in_known, na.rm=TRUE)

                            if(is.na(delta_mz_cut) | delta_mz_cut < 3 * delta_mz){
                                delta_mz_cut <- 3 * delta_mz
                            }

                            dif1 <- detected_peaks_rt_corrected - mean(RT_in_known_corrected, na.rm=TRUE)
                            dif2 <- detected_peaks_mz_corrected - mean(mz_in_known_corrected, na.rm=TRUE)
                            within_range <- ifelse(abs(dif1) <= delta_RT_cut & abs(dif2) <= delta_mz_cut, TRUE, FALSE)

                            if(any(within_range)){
                                detected_peaks_rt <- detected_peaks_rt[within_range]
                                detected_peaks_mz <- detected_peaks_mz[within_range]
                                detected_peaks_rt_corrected <- detected_peaks_rt_corrected[within_range]
                                detected_peaks_mz_corrected <- detected_peaks_mz_corrected[within_range]

                                # Now check that at the position of the
                                # potential peak no other peak is present in
                                # samples with known peak
                                count_other_peaks <- sum(peak_quant$Peak_1$num_peaks_with_background[i, ..kpi] - 1)
                                if(count_other_peaks > 0){
                                    other_peaks <- base::as.data.frame(matrix(ncol=6, nrow=count_other_peaks, 0))
                                    counter <- 0

                                    # Collect peak data for all additional peaks
                                    # in samples were peak was known
                                    for(o in kpi){
                                        cur_peak_count <- as.numeric(peak_quant$Peak_1$num_peaks_with_background[i, ..o] - 1)

                                        if(cur_peak_count > 0){
                                            for(p in 1:cur_peak_count){
                                                counter <- counter + 1
                                                temp <- c(as.numeric(peak_quant[[p + 1]]$Peak_mz_with_background[i, ..o] - mz_correction_factors[i, o]),
                                                          as.numeric(peak_quant[[p + 1]]$Peak_rt_with_background[i, ..o] - RT_correction_factors[i, o]),
                                                          as.numeric(peak_quant[[p + 1]]$feat_w_backgr_int[i, ..o]),
                                                          as.numeric(peak_quant[[p + 1]]$Icount_feat_w_bkgr_int[i, ..o]))
                                                # Add distance in RT and mz to
                                                # known peak
                                                temp[5:6] <- c(abs(temp[1] - as.numeric(peak_quant[[1]]$Peak_mz_with_background[i, ..o] - mz_correction_factors[i, o])),
                                                               abs(temp[2] - as.numeric(peak_quant[[1]]$Peak_rt_with_background[i, ..o] - RT_correction_factors[i, o])))
                                            }
                                        }
                                    }
                                    # Disregard other peaks below significance
                                    # threshold
                                    other_peaks <- other_peaks[which(other_peaks$V4 > peak_min_ion_count), ]
                                    # Disregard other peaks which are too close
                                    # to the expected peak
                                    # Cutoff RT: > 1.1* half peak width -> half RT extraction window
                                    # Cutoff mz: > 1.1* delta_mz -> half mz extraction window
                                    other_peaks <- other_peaks[which(other_peaks$V5 > delta_mz * 1.1 | other_peaks$V6 > ((features_select$RT_length[i] / 2) * 1.1)), ]
                                    # Now check that selected peaks are not close to other peaks in sampels with known peak
                                    dif1 <- detected_peaks_rt_corrected - other_peaks$V2
                                    dif2 <- detected_peaks_mz_corrected - other_peaks$V1
                                    overlap <- ifelse(any(abs(dif1) <= delta_rt / 2 & abs(dif2) <= delta_mz / 2, na.rm=TRUE), TRUE, FALSE)
                                }

                                else{
                                    overlap <- FALSE
                                }

                                # Not overlapping with any other peak
                                if(any(overlap == FALSE)){
                                    selected_peak <- 1

                                    data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant[[selected_peak]]$features_intensity[i, ..s]))
                                    data.table::set(Ioncount_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant[[selected_peak]]$Icount_feat_sample_mat[i, ..s]))
                                    data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant[[selected_peak]]$feat_w_backgr_int[i, ..s]))
                                    data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant[[selected_peak]]$Icount_feat_w_bkgr_int[i, ..s]))
                                    data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(selected_peak))
                                }

                                # The peak is overlapping with an other peak
                                # detected in other samples where correct peak
                                # was known
                                else{
                                    data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$features_intensity[i, ..s]))                                    data.table::set(Ioncount_sample,i = as.integer(i),j = as.integer(1),value = as.numeric(peak_quant$Standard$Icount_feat_sample_mat[i,..s]))
                                    data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$feat_w_backgr_int[i, ..s]))
                                    data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_w_bkgr_int[i, ..s]))
                                    data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(length(peak_quant)))
                                }
                            }

                            # None of the peaks is within the accepted range
                            # thus use standard window
                            else{
                                data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$features_intensity[i, ..s]))
                                data.table::set(Ioncount_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_sample_mat[i, ..s]))
                                data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$feat_w_backgr_int[i, ..s]))
                                data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_w_bkgr_int[i, ..s]))
                                data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(length(peak_quant)))
                            }
                        }

                        # No peak was valid
                        else{
                            data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$features_intensity[i, ..s]))
                            data.table::set(Ioncount_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_sample_mat[i, ..s]))
                            data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$feat_w_backgr_int[i, ..s]))
                            data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_w_bkgr_int[i, ..s]))
                            data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(length(peak_quant)))
                        }
                    }

                    # In no other sample the correct peak was known thus use the
                    # standard window
                    else{
                        data.table::set(features_intensity_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$features_intensity[i, ..s]))
                        data.table::set(Ioncount_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_sample_mat[i, ..s]))
                        data.table::set(feat_w_bkgrnd_inty_samp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$feat_w_backgr_int[i, ..s]))
                        data.table::set(Icount_w_bkgr_smp, i=as.integer(i), j=as.integer(1), value=as.numeric(peak_quant$Standard$Icount_feat_w_bkgr_int[i, ..s]))
                        data.table::set(peak_selected_sample, i=as.integer(i), j=as.integer(1), value=as.numeric(length(peak_quant)))
                    }
                }
            }
        }

        updatecounter <- updatecounter + 1
        if(updatecounter >= 10 & progress == TRUE){
            time_elapsed <- difftime(Sys.time(), start_time,units="secs")
            time_require <- (time_elapsed / (i / max)) * (1 - (i / max))
            td <- lubridate::seconds_to_period(time_require)
            time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                    lubridate::minute(td),
                                    round(lubridate::second(td), digits=0))
            updatecounter <- 0
            label <- base::paste(round(i / max * 100, 0), " % done (", i, "/",
                                 max, ", Time require: ", time_require, ")",
                                 sep="")
            tcltk::setTkProgressBar(pb, i, label=label)
        }
    }

    if(progress == TRUE){
        close(pb)
    }

    return(list(features_intensity_sample=features_intensity_sample,
                Ioncount_sample=Ioncount_sample,
                feat_w_bkgrnd_inty_samp=feat_w_bkgrnd_inty_samp,
                Icount_w_bkgr_smp=Icount_w_bkgr_smp,
                peak_selected_sample=peak_selected_sample))
}

#' Peak selection FDR algorithm for internal use
#' @param num_features num_features
#' @param features features
#' @param samples samples
#' @param peaks peaks
#' @param path_to_features path_to_features
#' @param peak_quant peak_quant
#' @param feat_w_backgr_int feat_w_backgr_int
#' @param peak_selected peak_selected
#' @param delta_mz delta_mz
#' @param delta_rt delta_rt
#' @param peak_min_ion_count peak_min_ion_count
#' @param num_peaks_store num_peaks_store
#' @param align_score_thr align_score_thr
#' @param n_cores n_cores
#' @param peak_decision peak_decision
#' @param Alignment_scoring Alignment_scoring
#' @param seed seed
#' @param plot plot
#' @import data.table
#' @import foreach
#' @import randomForest
#' @export
#' @details Peak selection FDR algorithm for internal use
Peak_selection_FDR <- function(num_features=500, features, samples, peaks,
                               path_to_features, peak_quant,
                               feat_w_backgr_int, peak_selected,
                               delta_mz, delta_rt, peak_min_ion_count,
                               num_peaks_store, align_score_thr,
                               n_cores, peak_decision, Alignment_scoring,
                               seed=110519, plot=TRUE){
    if(length(which(peaks$known == 1)) >= 10){
        if(!is.na(seed)){
            set.seed(seed)
        }

        # Select features for FDR estimation per sample
        c1 <- peaks$peak == peaks$selected
        c2 <- peaks$known == 1
        c3 <- !grepl("_i|_d", peaks$feature)
        temp_peaks <- peaks[which(c1 & c2 & c3), ]
        selection_feature <- sample(unique(temp_peaks$feature),
                                    length(unique(temp_peaks$feature)),
                                    replace=FALSE)
        selection_feature <- match(selection_feature, features$Feature_name)
        selection_feature_per_sample <- list()

        for(s in 1:length(samples)){
            temp_peaks_2 <- temp_peaks[which(temp_peaks$sample == samples[s]), ]
            userows <- match(features$Feature_name[selection_feature],
                             temp_peaks_2$feature)
            temp_peaks_2 <- temp_peaks_2[userows, ]
            sel <- which(!is.na(temp_peaks_2$known))[1:num_features]

            if(any(is.na(sel))){
                sel <- sel[which(!is.na(sel))]
            }

            selection_feature_per_sample[[s]] <- selection_feature[sel]
        }

        # Remove known peak location for selected features and samples and
        # replace with predicted corrections based on models
        features_select_FDR <- list()
        features$Observed_mz <- as.character(features$Observed_mz)
        features$Observed_RT <- as.character(features$Observed_RT)

        setwd(base::paste(path_to_features, "/Temporary_files", sep=""))
        e <- new.env()
        load("Feature_alignment_QC_data.RData", envir=e)
        qc <- e$QC_data
        median_feature_properties <- qc$mz_calibration_median_feature_properties
        mz_correction_models <- qc$mz_calibration_models
        RT_alignment_GAM_models <- qc$RT_calibration_GAM_models

        for(s in 1:length(samples)){
            temp <- data.table::copy(features)
            temp <- temp[selection_feature_per_sample[[s]], ]

            if(nrow(temp) > 10){
                # Remove observed RT and mz
                obs_rt <- (stringr::str_split(temp$Observed_RT, ";",
                                              simplify=TRUE))
                obs_rt[, s] <- "NA"
                obs_mz <- (stringr::str_split(temp$Observed_mz, ";",
                                              simplify=TRUE))
                obs_mz[, s] <- "NA"
                temp$Observed_mz <- apply(obs_mz, 1, base::paste, collapse=";")

                # Exchange mz calibration with prediction from model
                select_model <- which(names(mz_correction_models) == samples[s])
                features_select <- selection_feature_per_sample[[s]]
                loc <- match(features_select, median_feature_properties$Feature)
                temp_data <- median_feature_properties[loc, ]
                temp_data$Resolution <- 0
                usecols <- c("Retention.time", "m.z", "Charge", "Resolution")
                loc <- which(rowSums(is.na(temp_data[, usecols])) == 0)
                temp_data <- temp_data[loc, ]

                if(length(select_model) > 0){
                    if(nrow(temp_data) > 0){
                        mod <- mz_correction_models[[select_model]]
                        prediction <- stats::predict(mod, temp_data[, -1],
                                                     type="response")

                        if(any(is.na(prediction))){
                            prediction[which(is.na(prediction))] <- 0
                        }

                        aux <- base::paste("mz_calibration.", samples[s],
                                           sep="")
                        loc <- which(colnames(temp) == aux)
                        data.table::set(temp,
                                        as.integer(match(rownames(temp_data),
                                                         rownames(temp))),
                                        as.integer(loc), prediction)
                    }

                    else{
                        aux <- base::paste("mz_calibration.", samples[s],
                                           sep="")
                        loc <- which(colnames(temp) == aux)
                        data.table::set(temp,
                                        as.integer(match(rownames(temp_data),
                                                         rownames(temp))),
                                        as.integer(loc), 0)
                    }
                }

                else{
                    data.table::set(temp,
                                    as.integer(match(rownames(temp_data),
                                                     rownames(temp))),
                                    as.integer(which(grepl("mz_calibration.",
                                                           colnames(temp)))[s]),
                                    0)
                }

                # Exchange RT calibration with prediction from model
                select_model <- s
                if(!is.null(RT_alignment_GAM_models[[select_model]])){
                    mod <- RT_alignment_GAM_models[[select_model]]
                    prediction <- stats::predict(mod,
                                                 base::data.frame(x=temp$RT))
                    if(any(is.na(prediction))){
                        prediction[which(is.na(prediction))] <- 0
                    }

                    aux <- base::paste("RT_calibration.", samples[s], sep="")
                    loc <- which(colnames(temp) == aux)
                    data.table::set(temp, as.integer(1:nrow(temp)),
                                    as.integer(loc), prediction)
                }

                else{
                    aux <- base::paste("RT_calibration.", samples[s], sep="")
                    loc <- which(colnames(temp) == aux)
                    data.table::set(temp, as.integer(1:nrow(temp)),
                                    as.integer(loc), 0)
                }

                features_select_FDR[[s]] <- temp
            }

            else{
                features_select_FDR[[s]] <- temp
            }
        }

        # Prepare columns containing correction information
        indices_RT_correction <- which(grepl("RT_calibration",
                                             colnames(features)))
        indices_mz_correction <- which(grepl("mz_calibration",
                                             colnames(features)))
        cols <- colnames(features)[indices_RT_correction]
        ordering_indices <- match(base::gsub("-", ".", samples),
                                  base::gsub("RT_calibration\\.", "",  cols))
        indices_RT_correction <- indices_RT_correction[ordering_indices]
        indices_mz_correction <- indices_mz_correction[ordering_indices]

        # Reduce peak_quant information to only cover features selected for FDR
        # calculation

        # Prepare only required peak quant data
        peak_quant_reduced <- data.table::copy(peak_quant)
        all_selected_features <- unique(unlist(selection_feature_per_sample))

        for(p in 1:(num_peaks_store + 1)){
            for(p2 in 1:length(peak_quant_reduced[[p]])){
                pk <- peak_quant_reduced[[p]][[p2]]
                peak_quant_reduced[[p]][[p2]] <- pk[all_selected_features, ]
            }
        }

        cl <- snow::makeCluster(n_cores)
        doSNOW::registerDoSNOW(cl)

        start <- Sys.time()

        res <- foreach::foreach(s=1:length(samples)) %dopar% {
            selection_feature <- selection_feature_per_sample[[s]]

            if(length(selection_feature) > 0){
                max <- 1
                pct <- round(0 / max * 100, 0)
                pb <- tcltk::tkProgressBar(title="Evaluate peak selection FDR",
                                           label=base::paste(pct, "% done"),
                                           min=0, max=max, width=300)

                ft_sel_FDR_cur <- features_select_FDR[[s]]

                # Now perform peak decision for FDR_features and indicate if
                # same peak was selected as before
                len <- length(selection_feature)
                pk_dec_sme_wo_pk_knwn <- vector("logical", len)
                RT_correction_factors <- ft_sel_FDR_cur[, indices_RT_correction]
                mz_correction_factors <- ft_sel_FDR_cur[, indices_mz_correction]
                colnames(RT_correction_factors) <- samples
                rownames(RT_correction_factors) <- ft_sel_FDR_cur$Feature_name
                colnames(mz_correction_factors) <- samples
                rownames(mz_correction_factors) <- ft_sel_FDR_cur$Feature_name

                # Alias
                ft_w_bk_ity <- feat_w_backgr_int

                ft_w_bkgrd_inty_sel_FDR <- data.table::copy(ft_w_bk_ity)
                peak_selected_select_FDR <- data.table::copy(peak_selected)

                name <- rep(0L, length(selection_feature))
                temp_dummy_df <-  data.table::data.table(Name=name)
                colnames(temp_dummy_df) <- samples[s]
                close(pb)

                peak_quant_temp <- data.table::copy(peak_quant_reduced)

                for(p in 1:(num_peaks_store + 1)){
                    for(p2 in 1:length(peak_quant_temp[[p]])){
                        pk <- peak_quant_temp[[p]][[p2]]
                        loc <- match(selection_feature, all_selected_features)
                        peak_quant_temp[[p]][[p2]] <- pk[loc, ]
                    }
                }

                peak_quant_temp$Peak_1$correct_peak_w_background[, s] <- 0L

                usefeat <- ft_w_bkgrd_inty_sel_FDR[selection_feature, s,
                                                   drop=FALSE]
                usepeak <- peak_selected_select_FDR[selection_feature, s,
                                                    drop=FALSE]
                res_temp <- peak_decision(
                    features_select=ft_sel_FDR_cur,
                    peak_quant=peak_quant_temp,
                    samples=samples,
                    s=s,
                    RT_correction_factors=RT_correction_factors,
                    mz_correction_factors=mz_correction_factors,
                    features_intensity_sample=temp_dummy_df,
                    Ioncount_sample=temp_dummy_df,
                    feat_w_bkgrnd_inty_samp=usefeat,
                    Icount_w_bkgr_smp=temp_dummy_df,
                    peak_selected_sample=usepeak,
                    delta_mz=delta_mz,
                    delta_rt=delta_rt,
                    peak_min_ion_count=peak_min_ion_count,
                    progress=TRUE
                )

                c1 <- res_temp$peak_selected_sample[, 1]
                c2 <- peak_selected[selection_feature, s]
                pk_dec_sme_wo_pk_knwn <- c1 == c2

                usefeat <- res_temp$feat_w_bkgrnd_inty_samp
                data.table::set(ft_w_bkgrd_inty_sel_FDR,
                                as.integer(selection_feature),
                                as.integer(s), usefeat)
                data.table::set(peak_selected_select_FDR,
                                as.integer(selection_feature),
                                as.integer(s), res_temp$peak_selected_sample)

                samp <- base::paste(ft_sel_FDR_cur$Feature_name, "_Sample_", s,
                                    sep="")
                names(pk_dec_sme_wo_pk_knwn) <- samp

                # peak_decision_same_without_peak_known = pk_dec_sme_wo_pk_knwn
                return(list(pk_dec_sme_wo_pk_knwn=pk_dec_sme_wo_pk_knwn,
                            ft_w_bkgrd_inty_sel_FDR=ft_w_bkgrd_inty_sel_FDR,
                            peak_selected_select_FDR=peak_selected_select_FDR))
            }

            else{
                return(NULL)
            }
        }

        snow::stopCluster(cl)

        # Determine how many false selections show large intensity difference
        # and would not be removed by alignment scoring
        results_peak_selection_FDR_all <- list()
        max <- length(samples)
        prg <- base::paste(round(0 / max * 100, 0), "% done")
        pb <- tcltk::tkProgressBar(title="Evaluate peak selection FDR",
                                   label=prg, min=0, max=max, width=300)
        start_time <- Sys.time()
        updatecounter <- 0
        time_require <- 0

        for(s in 1:length(samples)){
            if(!is.null(res[[s]])){
                pk_dec_sme_wo_pk_knwn <- res[[s]]$pk_dec_sme_wo_pk_knwn
                peak_selected_select_FDR <- res[[s]]$peak_selected_select_FDR
                ft_w_bkgrd_inty_sel_FDR <- res[[s]]$ft_w_bkgrd_inty_sel_FDR
                selection_feature <- selection_feature_per_sample[[s]]

                if(length(which(pk_dec_sme_wo_pk_knwn == F)) > 0L){
                    selections_all <- 1:length(pk_dec_sme_wo_pk_knwn)

                    # Compare quantification differences for all tests
                    mat <- matrix(ncol=2, nrow=length(selections_all))
                    quant_res_4wrong <- base::as.data.frame(mat)
                    colnames(quant_res_4wrong) <- c("correct", "wrong")
                    ids <- features_select_FDR[[s]]$Feature_name[selections_all]
                    rownames(quant_res_4wrong) <- ids
                    corr <- as.numeric(quant_res_4wrong$correct)
                    quant_res_4wrong$correct <- corr
                    wrng <- as.numeric(quant_res_4wrong$wrong)
                    quant_res_4wrong$wrong <- wrng

                    for(i in selections_all){
                        ft <- selection_feature[i]
                        int_sel <- ft_w_bkgrd_inty_sel_FDR[ft, s]
                        ft_bck <- feat_w_backgr_int[ft, s]
                        data.table::set(quant_res_4wrong, as.integer(i),
                                        as.integer(1:2),
                                        list(base::log2(10^int_sel),
                                        base::log2(10^ft_bck)))
                    }
                    # Compare intensities of unbiased peak selection and correct
                    # peak
                    # Determine distribution of quantification deviations
                    # between correct and unbiased selected peak
                    if(plot == TRUE){
                        title <- base::paste(samples[s],
                                             "- Error in quantification")
                        graphics::plot(quant_res_4wrong$correct,
                                       quant_res_4wrong$wrong,
                                       xlab="True peak quantification, log2",
                                       ylab="Masked peak quantification, log2",
                                       main=title)
                        graphics::abline(a=0, b=1)
                        graphics::abline(a=1, b=1, lty=2, col="red")
                        graphics::abline(a=-1, b=1, lty=2, col="red")
                    }

                    # How many show deviation > 2 fold
                    dif <- quant_res_4wrong$wrong - quant_res_4wrong$correct
                    wrong_pk_inty_outlier <- abs(dif) > 1

                    # Finally visualize results
                    freq <- plyr::count(pk_dec_sme_wo_pk_knwn)

                    if(length(which(freq$x == FALSE)) == 0){
                        freq <- rbind(freq, base::data.frame(x=FALSE, freq=0))
                    }

                    if(length(which(freq$x == TRUE)) == 0){
                        freq <- rbind(freq, base::data.frame(x=TRUE, freq=0))
                    }

                    freq$rel <- freq$freq / sum(freq$freq) * 100

                    len <- length(which(wrong_pk_inty_outlier == TRUE))
                    plot_data <- c(freq$rel[which(freq$x == FALSE)],
                                   len / sum(freq$freq) * 100)
                    names(plot_data) <- c("Total",
                                          ">2-fold intensity difference")

                    results_peak_selection_FDR <- list(
                        pk_dec_sme_wo_pk_knwn=pk_dec_sme_wo_pk_knwn,
                        quant_res_4wrong=quant_res_4wrong,
                        wrong_pk_inty_outlier=wrong_pk_inty_outlier,
                        plot_data=plot_data
                    )
                }

                else{
                    freq <- plyr::count(pk_dec_sme_wo_pk_knwn)

                    if(length(which(freq$x == FALSE)) == 0){
                        freq <- rbind(freq, base::data.frame(x=FALSE, freq=0))
                    }

                    if(length(which(freq$x == TRUE)) == 0){
                        freq <- rbind(freq,base::data.frame(x=TRUE, freq=0))
                    }

                    freq$rel <- freq$freq / sum(freq$freq) * 100

                    plot_data <- c(freq$rel[which(freq$x == FALSE)],
                                   0 / sum(freq$freq) * 100)
                    names(plot_data) <- c("Total",
                                          ">2-fold intensity difference")

                    results_pk_selection_FDR <- list(
                        pk_dec_sme_wo_pk_knwn=pk_dec_sme_wo_pk_knwn,
                        quant_res_4wrong=NULL,
                        wrong_pk_inty_outlier=NULL,
                        plot_data=plot_data
                    )
                }

                results_peak_selection_FDR_all[[s]] <- results_pk_selection_FDR
            }

            else{
                results_peak_selection_FDR_all[[s]] <- NULL
            }

            updatecounter <- updatecounter + 1

            if(updatecounter >= 1){
                time_elapsed <- difftime(Sys.time(), start_time, units="secs")
                time_require <- (time_elapsed / (s / max)) * (1 - (s / max))
                td <- lubridate::seconds_to_period(time_require)
                time_require <- sprintf('%02d:%02d:%02d', td@hour,
                                        lubridate::minute(td),
                                        round(lubridate::second(td), digits=0))
                updatecounter <- 0
                lbl <- base::paste(round(s / max * 100, 0), " % done (", s,
                                   "/", max, ", Time require: ", time_require,
                                   ")", sep="")
                tcltk::setTkProgressBar(pb, s, label=lbl)
            }
        }

        close(pb)

        # Plot results over all samples
        Total_FDR <- NULL
        Large_Intensity_delta_FDR <- NULL

        for(s in 1:length(samples)){
            res <- results_peak_selection_FDR_all[[s]]

            # Require at least 100 identified peptides in a sample
            if(length(res$pk_dec_sme_wo_pk_knwn) >= 100){
                Total_FDR <- append(Total_FDR, res$plot_data[1])
                Large_Intensity_delta_FDR <- append(Large_Intensity_delta_FDR,
                                                    res$plot_data[2])
            }

            else{
                print(base::paste(samples[s], ": Not enough known peaks for ",
                                  "peak selection FDR estimation", sep=""))
                Total_FDR <- append(Total_FDR, NA)
                Large_Intensity_delta_FDR <- append(Large_Intensity_delta_FDR,
                                                    NA)
            }
        }

        names(Total_FDR) <- samples
        names(Large_Intensity_delta_FDR) <- samples

        if(plot == TRUE){
            ylim <- c(0, ifelse(max(Large_Intensity_delta_FDR, na.rm=TRUE) < 5,
                                5, max(Large_Intensity_delta_FDR, na.rm=TRUE)))
            title <- base::paste("Peak selection FDR - > 2-fold intensity ",
                                 "difference\nBased on n=", num_features,
                                 " random draws", sep="")
            p <- Barplots(Large_Intensity_delta_FDR, AvgLine=TRUE,
                          digits_average=1, Name=samples, xlab="",
                          ylab="FDR [%]", main=title, shownumbers=FALSE,
                          ylim=ylim)
            graphics::abline(h=5, lty=2, col="red")
        }
    }

    else{
        print(paste("Warning: The true peak for too few features is known. ",
                    "Skipping peak selection FDR estimation.", sep=""))
        results_peak_selection_FDR_all <- NA
        Total_FDR <- NA
        Large_Intensity_delta_FDR <- NA
    }

    return(list(results_peak_selection_FDR_all=results_peak_selection_FDR_all,
                Total_FDR=Total_FDR,
                Large_Intensity_delta_FDR=Large_Intensity_delta_FDR))
}
