
dropdownButton <- function(label="",
                           status=c("default", "primary", "success", "info",
                                    "warning", "danger"), ..., width=NULL){
    status <- match.arg(status)
    # Dropdown button content
    html_ul <- list(
        class="dropdown-menu",
        style=if(!is.null(width)){
            paste0("width: ", validateCssUnit(width), ";")
        },
        lapply(X=list(...), FUN=tags$li,
               style="margin-left: 10px; margin-right: 10px;")
    )
    # Dropdown button apparence
    html_button <- list(
        class=paste0("btn btn-", status," dropdown-toggle"),
        type="button",
        `data-toggle`="dropdown"
    )
    html_button <- c(html_button, list(label))
    html_button <- c(html_button, list(tags$span(class="caret")))

    # Final result
    tags$div(
        class="dropdown",
        do.call(tags$button, html_button),
        do.call(tags$ul, html_ul),
        tags$script("$('.dropdown-menu').click(function(e) {
                     e.stopPropagation();
                     });")
    )
}

# Panel for settings and run app
run_panel <- fluidPage(
    shinyjs::useShinyjs(),
    fluidRow(
        column(
            3,
            h3("Run"),
            actionButton("run_all", "Run", width=170),
            br(),
            radioButtons("massspecmode", h3("MassSpec-Mode"),
                         choices=list("Orbitrap"=1, "timsTOF"=2), selected=1),
            br()
        ),
        column(
            3,
            h3("Raw files"),
            actionButton("Raw_folder", "Choose directory", width=170),
            verbatimTextOutput("Raw_folder", placeholder=TRUE),
            helpText("Select the folder containing MS raw files")
        ),
        column(
            3,
            h3("MaxQ folder"),
            actionButton("MaxQ_output_folder", "Choose directory", width=170),
            verbatimTextOutput("MaxQ_output_folder", placeholder=TRUE),
            helpText(paste0("Select the folder containing MaxQuant output ",
                            "files e.g. the txt output folder"))
        ),
        column(
            3,
            h3("Results folder"),
            actionButton("IceR_output", "Choose directory", width=170),
            verbatimTextOutput("IceR_output", placeholder=TRUE),
            helpText("Select the folder where the IceR results should be saved")
        )
    ),
    fluidRow(
        column(
            3,
            textInput("Analysis_name", h3("Analysis name"),
                      value="IceR_analysis"),
            helpText("Specify how the analysis results should be named.")
        ),
        column(
            3,
            h3("Alignment settings")
        ),

        column(
            3,
            br(),
            br(),
            br(),
            h5("Minimal RT-Window [min]"),
            sliderInput("min_RT_window", "", min=0, max=10, value=1.0, step=0.1)
        ),
        column(
            3,
            br(),
            br(),
            br(),
            h5("Minimal m/z-Window [Da]"),
            sliderInput("min_mz_window", "", min=0, max=0.01, value=0.001,
                        step=0.0005)
        )
    ),
    fluidRow(
        column(
            3,
            h3("Load Settings"),
            actionButton("load_settings", "Load Settings", width=170)
        ),
        column(
            3,
            h3("Requantification settings")
        ),
        column(
            3,
            br(),
            br(),
            br(),
            h5("KDE resolution"),
            sliderInput("kde_resolution", "", min=10, max=200, value=50, step=5)
        ),
        column(
            3,
            br(),
            br(),
            br(),
            h5("Number of threads"),
            sliderInput("n_cores", "", min=1, max=50, value=8, step=1)
        )
    )
)

qc_panel <- fluidPage(
    fluidRow(
        column(
            2,
            h3("Visualize"),
            actionButton("vis_qc_alignment", "Alignment", width=100),
            br(),
            actionButton("vis_qc_quantification", "Quantification", width=100)
        ),
        column(
            5,
            plotOutput("plot_1", height="300px")
        ),
        column(
            5,
            plotOutput("plot_2", height="300px")
        )
    ),
    fluidRow(
        column(
            2
        ),
        column(
            5,
            plotOutput("plot_3", height="300px")
        ),
        column(
            5,
            plotOutput("plot_4", height="300px")
        )
    ),
    fluidRow(
        column(
            1
        ),
        column(
            5,
            plotOutput("plot_5", height="300px")
        ),
        column(
            6,
            plotOutput("plot_6", height="300px")
        )
    ),
    fluidRow(
        column(
            12,
            plotOutput("plot_7")
        )
    ),
    fluidRow(
        column(
            12,
            plotOutput("plot_8")
        )
    )
)


# Define UI ----
ui <- fluidPage(
    tags$style(type="text/css", ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"),
    titlePanel(strong("IceR")),
    tabsetPanel(type="tabs", tabPanel("Run", run_panel),
                tabPanel("QC", qc_panel))
)

run_all_processes <- function(settings_list){
    # Default parameters
    align_unknown <- FALSE
    only_unmodified <- FALSE
    RT_calibration <- TRUE
    mz_calibration <- TRUE
    use_isotopes <- TRUE
    intensity_correction <- TRUE
    peak_detection <- TRUE
    add_PMPs <- FALSE
    plot_peak_detection <- FALSE
    calc_protein_LFQ <- TRUE

    path_to_MaxQ_output <- ifelse(
        endsWith(settings_list$MaxQ_output_folder, "/"),
        substr(settings_list$MaxQ_output_folder, 1,
               nchar(settings_list$MaxQ_output_folder) - 1),
        settings_list$MaxQ_output_folder
    )



    path_to_output <- ifelse(endsWith(settings_list$IceR_output, "/"),
                             substr(settings_list$IceR_output, 1,
                                    nchar(settings_list$IceR_output) - 1),
                             settings_list$IceR_output)
    min_mz_window <- ifelse(settings_list$min_mz_window > 0,
                            settings_list$min_mz_window, NA)
    mz_window <- NA
    min_RT_window <- ifelse(settings_list$min_RT_window > 0,
                            settings_list$min_RT_window, NA)
    RT_window <- NA
    n_cores <- settings_list$n_cores
    MassSpec_mode <- ifelse(settings_list$MassSpec_settings == "1", "Orbitrap",
                            "TIMSToF")
    use_IM_data <- ifelse(settings_list$MassSpec_settings == "2", TRUE, FALSE)

    path_to_extracted_spectra <- NA
    path_to_mzXML <- NA

    if(MassSpec_mode == "TIMSToF"){
        path_to_extracted_spectra <- paste(settings_list$Raw_folder, "/",
                                           "all_ion_lists", sep="")
    }

    else{
        path_to_mzXML <- paste(settings_list$Raw_folder, "/", "mzXML", sep="")
    }

    max_steps <- 4

    if(use_isotopes == TRUE){
        max_steps <- max_steps + 1
    }

    if(add_PMPs == TRUE){
        max_steps <- max_steps + 1
    }

    cur_step <- 0

    withProgress(message='', min=0, max=max_steps, {
        print(paste0(Sys.time(), " Job started (IceR: ",
                     packageVersion("IceR"), ")"))
        print(R.version)

        # Save all parameters
        setProgress(cur_step, message="Prepare settings and parameters")

        Parameters <- as.data.frame(matrix(ncol=2, nrow=26))
        colnames(Parameters) <- c("Settings_name", "Setting")
        Parameters$Settings_name <- c("Path to raw files",
                                      "Path to MaxQ results", "Path to results",
                                      "Analysis name", "mz_window", "RT_window",
                                      "min_mz_window", "min_RT_window",
                                      "feature_mass_deviation_collapse",
                                      "only_unmodified_peptides",
                                      "align_unknown", "use_isotope_peaks",
                                      "peak_detection",
                                      "abundance_estimation_correction",
                                      "alignment_score_cut", "Quant_pVal_cut",
                                      "n_cores", "RT_correction",
                                      "mz_correction", "add_PMPs",
                                      "plot_peak_detection", "calc_protein_LFQ",
                                      "kde_resolution", "num_peaks_store",
                                      "MassSpec_mode", "use_IM_data")
        Parameters$Setting <- c(settings_list$Raw_folder, path_to_MaxQ_output,
                                path_to_output, settings_list$analysis_name,
                                mz_window, RT_window, min_mz_window,
                                min_RT_window, settings_list$collapse_mz,
                                only_unmodified, align_unknown, use_isotopes,
                                peak_detection, intensity_correction,
                                settings_list$Alignment_score_cut,
                                settings_list$Quant_pVal_cut, n_cores,
                                RT_calibration, mz_calibration, add_PMPs,
                                plot_peak_detection, calc_protein_LFQ,
                                settings_list$kde_resolution,
                                settings_list$num_peaks_store, MassSpec_mode,
                                use_IM_data)
        sample_list <- list.files(settings_list$Raw_folder)
        sample_list <- sample_list[which(grepl("\\.raw|\\.d", sample_list))]
        sample_list_process <- sample_list[which(!grepl("Library|Lib",
                                                        sample_list))]
        req <- ifelse(sample_list %in% sample_list_process, TRUE, FALSE)
        Raw_files <- data.frame(Raw_files=sample_list, Requantified=req)

        SaveExcel(list(Parameters=Parameters, Raw_files=Raw_files),
                  File=paste(path_to_output, "/Parameters.xlsx", sep=""))
        print(paste(Sys.time(), "Preparation finished"))

        # Convert raw files
        if(MassSpec_mode == "Orbitrap"){
            #Convert into mzXML
            cur_step <- cur_step + 1
            setProgress(cur_step, message="Convert raw files to mzXML")
            run_msconvert_raw_mzXML(settings_list$Raw_folder)
            print(paste(Sys.time(), "Raw file conversion finished"))
            ###Convert mzXML into usable RData files
            setProgress(cur_step, message="Prepare mzXML files")
            mzxml_to_list(path_to_mzXML, n_cores=n_cores)
            print(paste(Sys.time(), "Preparation of mzXMLs finished"))
        }

        else{
            cur_step <- cur_step + 1
            setProgress(cur_step, message="Extract MS1 spectra from raw files")
            convert_rawTIMS(path_to_raw=settings_list$Raw_folder)
            print(paste(Sys.time(), "Extraction finished"))
        }

        # Align features
        cur_step <- cur_step + 1
        setProgress(cur_step, message="Perform feature alignment")
        sample_list <- list.files(settings_list$Raw_folder)
        sample_list <- sample_list[which(grepl("\\.raw|\\.d", sample_list))]
        sample_list <- gsub("\\.raw|\\.d", "", sample_list)

        # Check that more than 1 sample is specified
        if(length(sample_list) < 2){
            stop(paste0("Less than 2 sample raw files are available in the ",
                        "Raw files folder."))
        }

        align_features(
            path_to_MaxQ_output=path_to_MaxQ_output,
            path_to_output=path_to_output,
            output_file_names_add=settings_list$analysis_name,
            mz_window=mz_window, RT_window=RT_window,
            min_mz_window=min_mz_window,
            min_RT_window=min_RT_window,
            min_num_ions_collapse=10,
            feature_mass_deviation_collapse=settings_list$collapse_mz,
            only_unmodified_peptides=only_unmodified,
            align_unknown=align_unknown,
            sample_list=sample_list,
            MassSpec_mode=MassSpec_mode,
            IM_window=NA,
            min_IM_window=0.002
        )

        print(paste(Sys.time(), "Feature alignmend finished"))

        # Add isotope features if required
        if(use_isotopes == TRUE){
            cur_step <- cur_step + 1
            setProgress(cur_step, message="Add +1 isotope features")
            ftfn <- paste("Features_aligned_merged_",
                          settings_list$analysis_name, ".txt", sep="")
            add_isotope_features(path_to_features=path_to_output,
                                 feature_table_file_name=ftfn,
                                 min_observations=0)
            print(paste(Sys.time(), "Addition of +1 isotope features finished"))
        }

        # Add potentially missed peptide features
        if(add_PMPs == TRUE){
            cur_step <- cur_step + 1
            setProgress(cur_step,
                        message="Add potentially missed peptide features")
            ftfn <- paste("Features_aligned_merged_",
                          settings_list$analysis_name, ".txt", sep="")
            add_missed_peptides(path_to_features=path_to_output,
                                feature_table_file_name=ftfn,
                                path_to_fasta=paste(path_to_MaxQ_output,
                                                    "/Search_DB.fasta", sep=""))
            print(paste(Sys.time(), "Addition of potentially missed peptide ",
                        "features finished"))
        }

        # Requantification of features
        cur_step <- cur_step + 1
        setProgress(cur_step,
                    message="Perform feature and protein quantification")
        ftfn <- paste("Features_aligned_merged_", settings_list$analysis_name,
                      ".txt", sep="")
        asc <- settings_list$Alignment_score_cut
        requantify_features(
            path_to_features=path_to_output,
            path_to_mzXML=path_to_mzXML,
            path_to_MaxQ_output=path_to_MaxQ_output,
            feature_table_file_name=ftfn,
            output_file_names_add=settings_list$analysis_name,
            RT_calibration=RT_calibration,
            mz_calibration=mz_calibration,
            n_cores=n_cores,
            abundance_estimation_correction=intensity_correction,
            Quant_pVal_cut=settings_list$Quant_pVal_cut,
            kde_resolution=settings_list$kde_resolution,
            num_peaks_store=settings_list$num_peaks_store,
            plot_peak_detection=plot_peak_detection,
            alignment_variability_score_cutoff=asc,
            alignment_scores_cutoff=settings_list$Alignment_score_cut,
            mono_iso_alignment_cutoff=settings_list$Alignment_score_cut,
            calc_protein_LFQ=calc_protein_LFQ,
            MassSpec_mode=MassSpec_mode,
            use_IM_data=use_IM_data,
            path_to_extracted_spectra=path_to_extracted_spectra
        )

        print(paste(Sys.time(), "Feature and protein quantification finished"))

        cur_step <- cur_step + 1
        setProgress(cur_step, message="Finishing")
    })

    print(paste(Sys.time(), "Job finished"))
}

# Server
server <- function(input, output, session){
    suppressWarnings(suppressMessages(library(IceR, quietly=TRUE)))
    suppressWarnings(suppressMessages(library(openxlsx, quietly=TRUE)))
    suppressWarnings(suppressMessages(library(shiny, quietly=TRUE)))
    suppressWarnings(suppressMessages(library(shinyFiles, quietly=TRUE)))
    suppressWarnings(suppressMessages(library(shinyjs, quietly=TRUE)))
    suppressWarnings(suppressMessages(library(stringr, quietly=TRUE)))
    suppressWarnings(suppressMessages(library(tcltk, quietly=TRUE)))
    suppressWarnings(suppressMessages(library(rChoiceDialogs, quietly=TRUE)))

    # Select file paths
    global_MaxQ_output_folder <- reactiveValues(datapath=getwd())
    global_Raw_folder <- reactiveValues(datapath=getwd())
    global_IceR_output <- reactiveValues(datapath=getwd())

    volumes <- getVolumes()

    MaxQ_output_folder <- reactive({
        return(parseDirPath(volumes, input$MaxQ_output_folder))
    })

    output$MaxQ_output_folder <- renderText({
    global_MaxQ_output_folder$datapath
    })

    observeEvent(input$MaxQ_output_folder, {
        if(Sys.info()["sysname"] != "Darwin"){
            selected_folder <- rchoose.dir(caption="Choose MaxQ output folder",
                                           default="~")
        }

        else{
            selected_folder <- tcltk::tk_choose.dir(
                caption="Choose MaxQ output folder",
                default="~"
            )
        }

        if(length(selected_folder) > 0){
            if(!is.na(selected_folder)){
                selected_folder <- gsub("\\\\", "/", selected_folder)
                global_MaxQ_output_folder$datapath <- selected_folder
            }
        }
    })

    # Choose Raw folder
    Raw_folder <- reactive({
        return(parseDirPath(volumes, input$Raw_folder))
    })

    output$Raw_folder <- renderText({
        global_Raw_folder$datapath
    })

    observeEvent(input$Raw_folder, {
        if(Sys.info()["sysname"] != "Darwin"){
            selected_folder <- rchoose.dir(
                caption="Choose folder containing raw files",
                default="~"
            )
        }

        else{
            selected_folder <- tcltk::tk_choose.dir(
                caption="Choose folder containing raw files",
                default="~"
            )
        }

        if(length(selected_folder) > 0){
            if(!is.na(selected_folder)){
                selected_folder <- gsub("\\\\", "/", selected_folder)
                global_Raw_folder$datapath <- paste0(selected_folder, "/")
            }
        }
    })

    # Choose IceR output folder
    IceR_output <- reactive({
        return(parseDirPath(volumes, input$IceR_output))
    })

    output$IceR_output <- renderText({
        global_IceR_output$datapath
    })

    observeEvent(input$IceR_output, {
        if(Sys.info()["sysname"] != "Darwin"){
            selected_folder <- rchoose.dir(caption="Choose final output folder",
                                           default="~")
        }

        else{
            selected_folder <- tcltk::tk_choose.dir(
                caption="Choose final output folder",
                default="~"
            )
        }

        if(length(selected_folder) > 0){
            if(!is.na(selected_folder)){
                selected_folder <- gsub("\\\\", "/", selected_folder)
                global_IceR_output$datapath <- selected_folder
            }
        }
    })

    # Combine all settings information into a list
    observeEvent(input$massspecmode, {
        #Orbitrap
        if(input$massspecmode == "1"){
            updateSliderInput(session, "min_mz_window", min=0, max=0.01,
                              value=0.001, step=0.0005)
        }

        else{
            updateSliderInput(session, "min_mz_window", min=0, max=0.01,
                              value=0.005, step=0.0005)
        }
    })

    # Buttons
    observeEvent(input$run_all, {
        settings_list <- list(
            MaxQ_output_folder=global_MaxQ_output_folder$datapath,
            Raw_folder=global_Raw_folder$datapath,
            IceR_output=global_IceR_output$datapath,
            min_RT_window=input$min_RT_window,
            RT_window=NA,
            min_mz_window=input$min_mz_window,
            mz_window=NA,
            collapse_mz=0.002,
            analysis_name=input$Analysis_name,
            Alignment_score_cut=0.05,
            Quant_pVal_cut=0.05,
            n_cores=input$n_cores,
            kde_resolution=input$kde_resolution,
            num_peaks_store=5,
            MassSpec_settings=input$massspecmode
        )
        run_all_processes(settings_list)
    })

    observeEvent(input$load_settings, {
        selected_parameters <- vector("character", 0)
        Filters_temp <- matrix(ncol=2, nrow=2)
        Filters_temp[1,] <- c("All files (*.*)", "*.*")
        Filters_temp[2,] <- c("Excel files (*.xlsx)", "*.xlsx")
        rownames(Filters_temp) <- c("All", "xlsx")

        rownames(Filters_temp)[nrow(Filters_temp)] <- "xlsx"

        if(Sys.info()["sysname"] != "Darwin"){
            selected_parameters <- rChoiceDialogs::rchoose.files(
                caption="Select IceR parameters file",
                multi=FALSE,
                filters=Filters_temp[c("All", "xlsx"), ],
                default="~"
            )
        }
        else{
            selected_parameters <- tcltk::tk_choose.files(
                caption="Select IceR parameters file",
                multi=FALSE,
                filters=Filters_temp[c("All", "xlsx"), ],
                default="~"
            )
        }


        if(length(selected_parameters) > 0){
            temp_settings <- read.xlsx(selected_parameters)
            temp_settings[is.na(temp_settings)] <- 0

            global_Raw_folder$datapath <- temp_settings$Setting[1]
            global_MaxQ_output_folder$datapath <- temp_settings$Setting[2]
            global_IceR_output$datapath <- temp_settings$Setting[3]

            updateTextInput(session, "Analysis_name",
                            value=temp_settings$Setting[4])
            updateSliderInput(session, "min_mz_window", min=0, max=0.01,
                              value=as.numeric(temp_settings$Setting[7]),
                              step=0.0005)
            updateSliderInput(session, "min_RT_window", min=0, max=10,
                              value=as.numeric(temp_settings$Setting[8]),
                              step=0.1)

            selection <- ifelse(temp_settings$Setting[10] == TRUE
                                & temp_settings$Setting[11] == TRUE, 1:2,
                                ifelse(temp_settings$Setting[10] == TRUE, 2,
                                       ifelse(temp_settings$Setting[11] == TRUE,
                                              1, 0)))

            selection <- NULL

            if(temp_settings$Setting[18] == TRUE){
                selection <- append(selection, 1)
            }

            if(temp_settings$Setting[19] == TRUE){
                selection <- append(selection, 2)
            }

            if(temp_settings$Setting[12] == TRUE){
                selection <- append(selection, 3)
            }

            if(temp_settings$Setting[14] == TRUE){
                selection <- append(selection, 4)
            }

            if(temp_settings$Setting[13] == TRUE){
                selection <- append(selection, 5)
            }

            if(temp_settings$Setting[20] == TRUE){
                selection <- append(selection, 6)
            }

            if(temp_settings$Setting[21] == TRUE){
                selection <- append(selection, 7)
            }

            if(temp_settings$Setting[22] == TRUE){
                selection <- append(selection, 8)
            }

            updateSliderInput(session, "n_cores", min=1, max=50,
                              value=as.numeric(temp_settings$Setting[17]),
                              step=1)
            updateSliderInput(session, "kde_resolution", min=10, max=200,
                              value=as.numeric(temp_settings$Setting[23]),
                              step=5)
            # Massspec mode
            selection_mode <- ifelse(temp_settings$Setting[25] == "TIMSToF"
                                     & temp_settings$Setting[26] == TRUE, 2, 1)
            updateRadioButtons(session, inputId="massspecmode",
                               choices=list("Orbitrap"=1, "timsTOF"=2),
                               selected=selection_mode)
        }
    })

    # Load QC data
    observeEvent(input$vis_qc_alignment, {
        path_to_QC_files <- paste(global_IceR_output$datapath,
                                  "/Temporary_files", sep="")
        load(paste(path_to_QC_files, "/Feature_alignment_QC_data.RData",
                   sep=""))

        # First row
        output$plot_1 <- renderPlot({
            visualize_MaxQ_RT_calibrations(QC_data)
        })

        output$plot_2 <- renderPlot({
            visualize_MaxQ_mz_calibrations(QC_data)
        })

        # Second row
        output$plot_3 <- renderPlot({
            aux <- QC_data$Alignment_deviations_overlap$sd_RT
            stats <- boxplot.stats(aux)$stats

            if(stats[5] > QC_data$Feature_alignment_windows$RT_window[2]){
                ylim <- c(0, stats[5])
            }

            else{
                ylim <- c(0, QC_data$Feature_alignment_windows$RT_window[2])
            }

            boxplot(aux, main="RT deviation between samples", outline=FALSE,
                    ylab="SD", ylim=ylim)
            abline(h=QC_data$Feature_alignment_windows$RT_window[2], lty=2,
                   col="red")
        })

        output$plot_4 <- renderPlot({
            aux <- QC_data$Alignment_deviations_overlap$sd_m.z
            stats <- boxplot.stats(aux)$stats

            if(stats[5] > QC_data$Feature_alignment_windows$mz_window[2]){
                ylim <- c(0, stats[5])
            }

            else{
                ylim <- c(0, QC_data$Feature_alignment_windows$mz_window[2])
            }

            boxplot(aux, main="m/z deviation between samples", outline=FALSE,
                    ylab="SD", ylim=ylim)
            abline(h=QC_data$Feature_alignment_windows$mz_window[2], lty=2,
                   col="red")
        })

        # Third row
        output$plot_5 <- renderPlot({
            aux <- QC_data$Alignment_deviations_overlap$overlap
            Barplots(plyr::count(aux)[, 2], main="Number of feature overlaps",
                     ylab="Count", Name=plyr::count(aux)[, 1], xlab="Overlaps",
                     AvgLine=FALSE, margins=c(6, 4, 4, 2))
        })

        output$plot_6 <- renderPlot({
            Visualize_mz_calibration_RF_models(QC_data)
        })

        # Fourth row
        output$plot_7 <- renderPlot({
            Visualize_RT_calibration(QC_data, path_to_QC_files)
        })

        # Fifth row
        output$plot_8 <- renderPlot({
            plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')
        })
    })

    observeEvent(input$vis_qc_quantification, {
        path_to_QC_files <- paste(global_IceR_output$datapath,
                                  "/Temporary_files", sep="")
        load(paste(path_to_QC_files, "/Feature_quantification_QC_data.RData",
                   sep=""))

        # First row
        output$plot_1 <- renderPlot({
            visualize_decoy_parameters_1(QC_data)
        })

        output$plot_2 <- renderPlot({
            visualize_decoy_parameters_2(QC_data)
        })

        # Second row
        output$plot_3 <- renderPlot({
            visualize_decoy_mean_intensities(QC_data)
        })

        output$plot_4 <- renderPlot({
            visualize_target_parameters(QC_data)
        })

        # Third row
        output$plot_5 <- renderPlot({
            visualize_peak_selection_FDR(QC_data)
        })

        output$plot_6 <- renderPlot({
            visualize_quantification_correction(QC_data)
        })

        # Fourth row
        output$plot_7 <- renderPlot({
           isualize_quant_significance(QC_data, 0.05)
        })

        # Fifth row
        output$plot_8 <- renderPlot({
            visualize_S2B(QC_data)
        })
    })
}

# Visualization functions
visualize_MaxQ_RT_calibrations <- function(QC_data){
    p <- boxplot(QC_data$MaxQ_calibrations$Retention.time.calibration,
                 main="Detected RT calibrations by MaxQ", ylab="RT correction",
                 outline=FALSE)
    text(1.3, p$stats[2], round(p$stats[2], digits=2))
    text(1.3, p$stats[3], round(p$stats[3], digits=2))
    text(1.3, p$stats[4], round(p$stats[4], digits=2))
    range <- p$stats[4] - p$stats[2]
}

visualize_MaxQ_mz_calibrations <- function(QC_data){
    p <- boxplot(QC_data$MaxQ_calibrations$Uncalibrated...Calibrated.m.z..Da.,
                 main="Detected m/z calibrations by MaxQ",
                 ylab="m/z correction", outline=FALSE)
    text(1.3, p$stats[2], round(p$stats[2], digits=4))
    text(1.3, p$stats[3], round(p$stats[3], digits=4))
    text(1.3, p$stats[4], round(p$stats[4], digits=4))
    range <- p$stats[4] - p$stats[2]
}

Visualize_mz_calibration_RF_models <- function(QC_data){
    rsq <- list()

    for(i in 1:length(names(QC_data$mz_calibration_models))){
        loc <- names(QC_data$mz_calibration_models)[i]
        rsq[[i]] <- QC_data$mz_calibration_models[[loc]]$rsq
    }

    rsq <- simplify2array(rsq)
    ylim <- c(min(rsq), max(rsq))
    colnames(rsq) <- names(QC_data$mz_calibration_models)

    for(i in 1:length(names(QC_data$mz_calibration_models))){
        if(i == 1){
            plot(rsq[, i], type="l", main="RF mz-calibration model - Rsq",
                 ylab="Rsq", xlab="Iterations", ylim=ylim)
            text(nrow(rsq), rsq[nrow(rsq), i], colnames(rsq)[i], pos=2)
        }

        else{
            lines(rsq[, i])
            text(nrow(rsq), rsq[nrow(rsq), i], colnames(rsq)[i], pos=2)
        }
    }
}

Visualize_RT_calibration <- function(QC_data, path_to_QC_files){
    library(matrixStats)

    e1 <- new.env()
    load(paste(path_to_QC_files,
               "/Quantification_raw_results_with_scores_filtered.RData",
               sep=""), e1)
    features <- e1$temp_results$features

    x_pred <- seq(min(features$RT, na.rm=TRUE), max(features$RT, na.rm=TRUE),
                  length.out=nrow(features))
    stats <- boxplot.stats(x_pred)
    x_pred <- x_pred[x_pred > stats$stats[2] & x_pred < stats$stats[4]]
    set.seed(1)
    x_pred <- sort(sample(x_pred, 1000))
    ncol <- length(QC_data$RT_calibration_GAM_models)
    deviation_from_0 <- as.data.frame(matrix(ncol=ncol, nrow=length(x_pred)))

    for(i in 1:length(QC_data$RT_calibration_GAM_models)){
        deviation_from_0[, i] <- predict(QC_data$RT_calibration_GAM_models[[i]],
                                         data.frame(x=x_pred))
    }

    colnames(deviation_from_0) <- names(QC_data$mz_calibration_models)
    median_devs <- colMedians(as.matrix(deviation_from_0), na.rm=TRUE)
    sd_devs <- colSds(as.matrix(deviation_from_0), na.rm=TRUE)
    ylim <- c(min(median_devs - sd_devs) - 2, max(median_devs + sd_devs) + 2)
    par(mar=c(8, 4, 4, 4))
    boxplot(deviation_from_0, outline=FALSE, ylab="RT-deviation [min]", las=2,
            ylim=ylim, main="RT-deviation from feature per sample")
    abline(h=0)
}

visualize_decoy_mean_intensities <- function(QC_data){
    library(mgcv)
    aux <- QC_data$Decoy_feature_parameters
    RT <- as.numeric(aux$Decoy_mean_intensity_per_RT$RT)
    x_pred <- seq(min(RT, na.rm=TRUE), max(RT, na.rm=TRUE),
                  length.out=length(RT))
    stats <- boxplot.stats(x_pred)
    x_pred <- x_pred[x_pred > stats$stats[2] & x_pred < stats$stats[4]]
    set.seed(1)
    x_pred <- sort(sample(x_pred, 1000))
    y_pred_all <- matrix(ncol=length(aux$GAM_model_per_sample), nrow=1000)

    for(i in 1:length(aux$GAM_model_per_sample)){
        fit_gam_mean <- aux$GAM_model_per_sample[[i]]
        y_pred_all[, i] <- predict(fit_gam_mean, data.frame(RT=x_pred))
    }

    ylim <- c(min(y_pred_all, na.rm=TRUE) * 0.9,
              max(y_pred_all, na.rm=TRUE) * 1.1)

    for(i in 1:length(aux$GAM_model_per_sample)){
        if(i == 1){
            plot(x_pred, y_pred_all[, i], type="l", ylim=ylim,
                 ylab="Intensity, log2", xlab="RT [min]",
                 main="Decoy feature intensity GAMs")
            text(x_pred[1000], y_pred_all[1000, i],
                 names(aux$GAM_model_per_sample)[i], pos=2, cex=0.7)
        }

        else{
            lines(x_pred, y_pred_all[, i])
            text(x_pred[1000], y_pred_all[1000, i],
                 names(aux$GAM_model_per_sample)[i], pos=2, cex=0.7)
        }
    }
}

visualize_decoy_parameters_1 <- function(QC_data){
    par(mfrow=c(1, 2))
    aux <- QC_data$Decoy_feature_parameters
    boxplot(as.numeric(as.matrix(aux$decoy_ioncount)), outline=FALSE,
            ylab="Number of ions", main="Decoy - Ions")
    boxplot(as.numeric(as.matrix(aux$decoy_intensities)), outline=FALSE,
            ylab="Summed intensity, log2", main="Decoy - Sum int.")
}

visualize_decoy_parameters_2 <- function(QC_data){
    par(mfrow=c(1, 2))
    aux <- QC_data$Decoy_feature_parameters
    boxplot(as.numeric(as.matrix(aux$decoy_mean_intensity)), outline=FALSE,
            ylab="Mean intensity, log2", main="Decoy - Mean int.")
    boxplot(as.numeric(as.matrix(aux$decoy_sd_intensity)), outline=FALSE,
            ylab="SD of intensity, log2", main="Decoy - SD int.")
}

visualize_target_parameters <- function(QC_data){
    par(mfrow=c(1, 2))
    boxplot(as.numeric(QC_data$Target_ion_counts), outline=FALSE,
            main="Target - Ions", ylab="Number of ions")
    aux <- QC_data$Signal_to_background_target_decoy
    loc <- which(aux$Target_Decoy == "target")
    boxplot(as.numeric(as.matrix(aux$S2B[loc, ])), outline=FALSE,
            main="Target - S2B", ylab="Signal/Background, log2")
    abline(h=0, lty=2, col="red")
}

visualize_peak_selection_FDR <- function(QC_data){
    aux <- QC_data$FDR_peak_selection$Large_Intensity_delta_FDR
    ylim <- c(0, ifelse(max(aux, na.rm=TRUE) < 5, 5, max(aux, na.rm=TRUE)))
    p <- Barplots(aux, AvgLine=TRUE, digits_average=1, Name=names(aux), xlab="",
                  ylab="FDR [%]", main="Peak selection FDR", shownumbers=FALSE,
                  ylim=ylim)
    abline(h=5, lty=2, col="red")
}

visualize_quantification_correction <- function(QC_data){
    dat <- QC_data$Abundance_correction$correction_data
    smoothScatter(dat[, 1],dat[, 2], ylab="Requant, log2", xlab="MaxQ, log2",
                  main="Correlation of MaxQ vs Requant on peptide level")
    abline(a=0, b=1)
    temp <- dat[c(order(dat[, 2], decreasing=TRUE)[1:(0.1 * nrow(dat))],
                  order(dat[, 2], decreasing=FALSE)[1:(0.1 * nrow(dat))]), ]
    y <- temp[, 2]
    x <- temp[, 1]
    fit <- lm(y~x)
    abline(fit, lty=2, col="red")
    m <- summary(fit)$coefficients[2, 1]
    Rsq <- summary(fit)$r.squared
    posx <- par("usr")[1] + (par("usr")[2] - par("usr")[1]) * 0.15
    posy <- par("usr")[4] - (par("usr")[4] - par("usr")[3]) * 0.1
    text(posx, posy, labels=paste("Rsq =", round(Rsq, digits=2), "\nslope =",
         round(m, digits=2)))
}

visualize_quant_significance <- function(QC_data, Quant_pVal_cut){
    pval_signal_with_background_quant <- QC_data$Quant_pval
    par(mar=c(10, 4, 4, 4))
    boxplot(-log10(pval_signal_with_background_quant), outline=FALSE,
            main="Quantification pVals per feature", ylab="pValue, -log10",
            names=colnames(pval_signal_with_background_quant), las=2)
    abline(h=-log10(Quant_pVal_cut), lty=2, col="red")
}

visualize_S2B <- function(QC_data){
    aux <- QC_data$Signal_to_background_target_decoy
    S2B <- aux$S2B[which(aux$Target_Decoy == "target"), ]
    par(mar=c(10, 4, 4, 4))
    boxplot(S2B, outline=FALSE,
            main="Signal to background ratio per feature quantification",
            ylab="Signal/Background, log2", names=colnames(S2B), las=2)
    abline(h=0, lty=2, col="red")

    for(co in 1:ncol(S2B)){
        med <- median(S2B[, co], na.rm=TRUE)
        text(co,med + ((par("usr")[4] - par("usr")[3]) * 0.05),
             round(med, digits=1))
    }
}

# Run the application
shinyApp(ui=ui, server=server)
