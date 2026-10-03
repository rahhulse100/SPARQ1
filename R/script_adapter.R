# Script-backed workflow adapter.

sparq_assess_script <-
function(data, script, input_object, result_object,
    output_type = c("scalar", "ranked", "partition", "graph"),
    working_directory = c("script", "temporary"), script_args = list(),
    ...) {

    if (!is.data.frame(data) || !nrow(data)) {
        stop("data must be a nonempty data.frame.", call. = FALSE)
    }
    if (!is.character(script) || length(script) != 1L || is.na(script) ||
        !nzchar(script) || !file.exists(script)) {
        stop("script must be one existing .R file.", call. = FALSE)
    }
    if (!is.character(input_object) || length(input_object) != 1L ||
        is.na(input_object) || !nzchar(input_object)) {
        stop("input_object must name the object read by the script.", call. = FALSE)
    }
    if (!is.character(result_object) || length(result_object) != 1L ||
        is.na(result_object) || !nzchar(result_object)) {
        stop("result_object must name the object created by the script.", call. = FALSE)
    }
    if (!is.list(script_args)) {
        stop("script_args must be a list.", call. = FALSE)
    }

    output_type <- match.arg(output_type)
    working_directory <- match.arg(working_directory)
    script_path <- normalizePath(script, mustWork = TRUE)

    run_script <- function(perturbed_data) {
        run_directory <- if (working_directory == "script") {
            dirname(script_path)
        } else {
            tempfile("sparq_script_")
        }

        if (working_directory == "temporary") {
            dir.create(run_directory, recursive = TRUE, showWarnings = FALSE)
            on.exit(unlink(run_directory, recursive = TRUE, force = TRUE), add = TRUE)
        }

        previous_directory <- getwd()
        on.exit(setwd(previous_directory), add = TRUE)
        setwd(run_directory)

        script_environment <- new.env(parent = .GlobalEnv)
        assign(input_object, perturbed_data, envir = script_environment)
        assign("sparq_script_args", script_args, envir = script_environment)

        sys.source(script_path, envir = script_environment)

        if (!exists(result_object, envir = script_environment, inherits = FALSE)) {
            stop(
                paste0(
                    "The script did not create `", result_object, "`. ",
                    "Assign its final analysis result to that object."
                ),
                call. = FALSE
            )
        }

        get(result_object, envir = script_environment, inherits = FALSE)
    }

    assessment <- sparq_assess_workflow(
        data = data,
        analysis_function = run_script,
        output_type = output_type,
        ...
    )

    assessment$script <- list(
        path = script_path,
        input_object = input_object,
        result_object = result_object,
        working_directory = working_directory,
        script_args = script_args
    )
    class(assessment) <- unique(c("sparq_script_assessment", class(assessment)))
    assessment
}
