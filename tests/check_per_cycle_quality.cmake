if(NOT EXISTS "${REPORT}")
    message(FATAL_ERROR
        "Per-cycle quality report not found: ${REPORT}")
endif()

file(STRINGS "${REPORT}" ROWS)

list(GET ROWS 0 HEADER)

if(NOT HEADER STREQUAL
   "cycle\tmean_quality\tlower_quartile\tmedian")
    message(FATAL_ERROR
        "Per-cycle quality report has an unexpected header: ${HEADER}")
endif()

list(LENGTH ROWS ROW_COUNT)

if(ROW_COUNT LESS 2)
    message(FATAL_ERROR
        "Per-cycle quality report contains no data rows")
endif()

set(PREVIOUS_END 0)

math(EXPR LAST_INDEX "${ROW_COUNT} - 1")

foreach(INDEX RANGE 1 ${LAST_INDEX})

    list(GET ROWS ${INDEX} ROW)

    string(REPLACE "\t" ";" FIELDS "${ROW}")
    list(LENGTH FIELDS FIELD_COUNT)

    if(NOT FIELD_COUNT EQUAL 4)
        message(FATAL_ERROR
            "Invalid number of columns in row: ${ROW}")
    endif()

    list(GET FIELDS 0 CYCLE)
    list(GET FIELDS 1 MEAN)
    list(GET FIELDS 2 LOWER_QUARTILE)
    list(GET FIELDS 3 MEDIAN)

    # -------------------------
    # Проверка cycle / BaseGroup
    # -------------------------

    if(CYCLE MATCHES "^([0-9]+)-([0-9]+)$")

        set(START "${CMAKE_MATCH_1}")
        set(END "${CMAKE_MATCH_2}")

    elseif(CYCLE MATCHES "^[0-9]+$")

        set(START "${CYCLE}")
        set(END "${CYCLE}")

    else()

        message(FATAL_ERROR
            "Invalid cycle/BaseGroup label: ${CYCLE}")

    endif()

    if(START LESS 1)
        message(FATAL_ERROR
            "Cycle range must start at >= 1: ${CYCLE}")
    endif()

    if(END LESS START)
        message(FATAL_ERROR
            "Invalid cycle range: ${CYCLE}")
    endif()

    # Группы должны идти строго вперёд.
    #
    # Допускаем 1, 2, 3, 14-15 и т.д.,
    # но не перекрывающиеся или повторяющиеся группы.

    if(START LESS_EQUAL PREVIOUS_END)
        message(FATAL_ERROR
            "Cycle groups overlap or are not ordered: ${CYCLE}")
    endif()

    set(PREVIOUS_END "${END}")

    # -------------------------
    # Проверка mean_quality
    # -------------------------

    if(NOT MEAN MATCHES "^[0-9]+(\\.[0-9]+)?$")
        message(FATAL_ERROR
            "Invalid mean_quality in row: ${ROW}")
    endif()

    # -------------------------
    # Проверка quartile
    # -------------------------

    if(NOT LOWER_QUARTILE MATCHES
       "^([0-9]+(\\.[0-9]+)?|nan)$")

        message(FATAL_ERROR
            "Invalid lower_quartile in row: ${ROW}")

    endif()

    # -------------------------
    # Проверка median
    # -------------------------

    if(NOT MEDIAN MATCHES
       "^([0-9]+(\\.[0-9]+)?|nan)$")

        message(FATAL_ERROR
            "Invalid median in row: ${ROW}")

    endif()

endforeach()