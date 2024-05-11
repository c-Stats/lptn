# nolint start

list.of.packages <- c("quantmod", "stringr", "data.table", "dplyr", "magrittr", "ggplot2", "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if (!"htmltab" %in% installed.packages()[,"Package"]){
    devtools::install_version("htmltab", version = "0.8.2", repos = "http://cran.us.r-project.org")
}

# Load previously installed packages
list.of.packages <- c(list.of.packages, "htmltab")
for (package in list.of.packages) {
    if (package %in% installed.packages()) {
        library(package, character.only = TRUE)
    }
}


n_days_histo = 10*365
path <- paste("./Data")
current_date <- Sys.Date()

#Create directory if it doesn't exist
dir.create(path, showWarnings = FALSE)

#Check if the SPY stock list is in the directory
SPY_stocks_path <- paste(path, "SPY_stocks.csv", sep = "/")
SPY_exists <- file.exists(SPY_stocks_path)

old_SPY_stocks_path <- paste(path, "old_SPY_stocks.csv", sep = "/")
old_SPY_exists <- file.exists(old_SPY_stocks_path)	

SPY_table_path <- paste(path, "SPY_table.csv", sep = "/")

if(!SPY_exists | !old_SPY_exists){

    message("SPY_stocks.csv not found. Downloading file from Wikipedia.")
    SPY_stocks <- htmltab::htmltab("https://en.wikipedia.org/wiki/List_of_S%26P_500_companies")
    write.csv(SPY_stocks, SPY_stocks_path, row.names = FALSE)
    SPY_stocks <- data.table::as.data.table(SPY_stocks)

    old_SPY_stocks <- htmltab::htmltab("https://en.wikipedia.org/wiki/List_of_S%26P_500_companies", 2)
    write.csv(old_SPY_stocks, old_SPY_stocks_path, row.names = FALSE)
    old_SPY_stocks <- data.table::as.data.table(old_SPY_stocks)


} else {

    SPY_stocks <- data.table::fread(SPY_stocks_path)
    old_SPY_stocks <- data.table::fread(old_SPY_stocks_path)

}

SPY_tickers <- unique(c(SPY_stocks$Symbol, old_SPY_stocks$`Removed >> Ticker`, old_SPY_stocks$`Added >> Ticker`))
SPY_tickers <- SPY_tickers[which(!is.na(SPY_tickers))]

#Check if the history is available 
#Create directory if needed
dir.create(file.path(path, "Quotes"), showWarnings = FALSE)

SPY_quotes_path <- paste(path, "Quotes", sep = "/")
SPY_exists <- file.exists(paste(SPY_quotes_path, "SPY_close.csv", sep = "/"))

if(!SPY_exists){

    from <- as.Date(current_date) - n_days_histo
    to <- current_date

    #Retrieve history
    message("No previous history found. Downloading stock quotes.")
    quantmod::getSymbols(SPY_tickers, 
                            from = from,
                            to = to,
                            warnings = FALSE, auto.assign = TRUE)


    #Flush non-loaded files
    SPY_tickers <- SPY_tickers[which(sapply(SPY_tickers, exists))]
    SPY_stocks <- SPY_stocks[Symbol %in% SPY_tickers]
    old_SPY_stocks <- old_SPY_stocks[eval(`Removed >> Ticker`) %in% SPY_tickers | eval(`Added >> Ticker`) %in% SPY_tickers]

    SPY_table <- data.table::copy(SPY_stocks[, c("Symbol", "Security", "Date added"), with = FALSE]) %>%
                                    .[, "Date removed" := NA]

    temp1 <- data.table::copy(old_SPY_stocks[, c("Removed >> Ticker", "Removed >> Security"), with = FALSE]) %>%
                                    .[, "Date removed" := as.Date(sapply(old_SPY_stocks$Date, function(x){as.Date(x, "%B %d, %Y")}), origin = "1970-01-01" )] 

    temp2 <- data.table::copy(old_SPY_stocks[, c("Added >> Ticker", "Added >> Security"), with = FALSE]) %>%
                                    .[, "Date added" := as.Date(sapply(old_SPY_stocks$Date, function(x){as.Date(x, "%B %d, %Y")}), origin = "1970-01-01")] 

    names(SPY_table) <- c("Symbol", "Security", "Added", "Removed")
    names(temp1) <- c("Symbol", "Security", "Removed")
    names(temp2) <- c("Symbol", "Security", "Added")

    setkey(temp1, "Symbol")
    setkey(temp2, "Symbol")	

    temp1[temp2, Added := i.Added]	
    temp2[temp1, Removed := i.Removed]	

    temp1 <- temp1[, c("Symbol", "Security", "Added", "Removed"), with = FALSE]	%>%
                    .[, Added := as.Date(Added)] %>%
                    .[, Removed := as.Date(Removed)]

    temp2 <- temp2[, c("Symbol", "Security", "Added", "Removed"), with = FALSE]	%>%
                    .[, Added := as.Date(Added)] %>%
                    .[, Removed := as.Date(Removed)]

    SPY_table <- SPY_table %>%
                    .[, Added := as.Date(Added)] %>%
                    .[, Removed := as.Date(Removed)]					

    SPY_table <- data.table::copy(dplyr::bind_rows(list(a = SPY_table, b = temp1, c = temp2))) %>%
                    .[, lapply(.SD, function(x){x[1]}), by = "Symbol", .SDcols = c("Security", "Added", "Removed")]


    #Group individual files into a proper list
    stocks <- lapply(SPY_tickers, get)
    names(stocks) <- SPY_tickers

    #Filter 
    classes <- unlist(lapply(stocks, function(x){class(x)[1]}))
    keep <- which(classes == "xts")
    SPY_tickers <- SPY_tickers[keep]
    stocks <- stocks[keep]
    SPY_stocks <- SPY_stocks[Symbol %in% SPY_tickers]

    #Discard files with missing dates
    getmode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }

    n_rows <- unlist(lapply(stocks, nrow))
    keep <- which(n_rows == getmode(n_rows))

    stocks <- lapply(stocks, as.data.table)

    template <- data.table::copy(stocks[[keep[1]]])
    to_na <- names(template)[-1]
    template[, eval(to_na) := rep(NA, nrow(template))]

    to_fix <- which(n_rows != getmode(n_rows))
    for(i in to_fix){

        dates_missing <- template[!(index %in% stocks[[i]]$index)]
        names(dates_missing) <- names(stocks[[i]])
        stocks[[i]] <- data.table::copy(dplyr::bind_rows(list(a = dates_missing, b = stocks[[i]])))
        setorder(stocks[[i]], "index")

    }

    n_rows <- unlist(lapply(stocks, nrow))
    keep <- which(n_rows == getmode(n_rows))

    stocks <- stocks[keep]
    SPY_tickers <- names(stocks)
    SPY_table <- SPY_table[Symbol %in% SPY_tickers]

    #Build Open, High, Low, Close, Volume, Adjusted frames
    stocks[-1] <- lapply(stocks[-1], function(x){if("index" %in% names(x)){x[, index := NULL]}})
    stocks <- data.table::copy(dplyr::bind_cols(stocks))

    varnames <- c("Open", "High", "Low", "Close", "Volume", "Adjusted")
    column_names <- names(stocks)

    quotes <- lapply(as.list(varnames), function(x){stocks[, c(1, which(grepl(x, column_names))), with = FALSE]})
    names(quotes) <- varnames
    for(i in 1:length(quotes)){

        names(quotes[[i]]) <- c("Date", SPY_tickers)
        quotes[[i]][, Date := as.Date(Date)]

        filename <- paste("SPY_", stringr::str_to_lower(names(quotes)[i]), ".csv", sep = "")
        path_save <- paste(SPY_quotes_path, filename, sep = "/")
        write.csv(quotes[[i]], path_save, row.names = FALSE)

    }

    SPY_table_path <- paste(path, "SPY_table.csv", sep = "/")
    write.csv(SPY_table, SPY_table_path, row.names = FALSE)

    message(paste("Successfully queried", nrow(SPY_stocks), "stocks from the SPDR S&P 500 ETF Trust."))	
    message(paste("Data saved at", SPY_quotes_path))

} else {

    varnames <- c("Open", "High", "Low", "Close", "Volume", "Adjusted")
    quotes <- lapply(as.list(varnames), function(x){

            filename <- paste("SPY_", stringr::str_to_lower(x), ".csv", sep = "")
            path_read <- paste(SPY_quotes_path, filename, sep = "/")
            data.table::fread(path_read)

        })

    names(quotes) <- varnames

    from <- max(quotes$Close$Date) 
    to <- current_date

    day_of_week <- weekdays(to)
    if(day_of_week == "Saturday"){

        to <- to - 1
        
    } else if(day_of_week == "Sunday"){

        to <- to - 2

    }

    #Update data if needed
    if(to > from){

        SPY_tickers <- unique(unlist(lapply(quotes, function(x){names(x)[-1]})))

        quantmod::getSymbols(SPY_tickers, 
                                from = from,
                                to = to + 1,
                                warnings = FALSE, auto.assign = TRUE)


        #Check if everything was loaded
        message("Pausing for 5 seconds...")
        Sys.sleep(5)
        was_loaded <- sapply(SPY_tickers, exists)
        #Re-try loading
        not_loaded <- which(!was_loaded)
        if(any(not_loaded)){

            n_not_loaded <- length(not_loaded)
            message(paste("Failed to load", n_not_loaded, "stock(s). Attempting to load them a second time."))

            quantmod::getSymbols(SPY_tickers[not_loaded], 
                                    from = from,
                                    to = to,
                                    warnings = FALSE, auto.assign = TRUE)

        }

        was_loaded <- sapply(SPY_tickers, exists)
        SPY_tickers <- SPY_tickers[was_loaded]

        #Group individual files into a proper list
        stocks <- lapply(SPY_tickers, get)
        names(stocks) <- SPY_tickers

        #Filter 
        classes <- unlist(lapply(stocks, function(x){class(x)[1]}))
        keep <- which(classes == "xts")
        stocks <- stocks[keep]
        SPY_tickers <- SPY_tickers[keep]

        #Discard files with missing dates
        getmode <- function(v) {
            uniqv <- unique(v)
            uniqv[which.max(tabulate(match(v, uniqv)))]
        }

        n_rows <- unlist(lapply(stocks, nrow))
        keep <- which(n_rows == getmode(n_rows))

        stocks <- lapply(stocks, as.data.table)

        template <- data.table::copy(stocks[[keep[1]]])
        to_na <- names(template)[-1]
        template[, eval(to_na) := rep(NA, nrow(template))]

        to_fix <- which(n_rows != getmode(n_rows))
        for(i in to_fix){

            dates_missing <- template[!(index %in% stocks[[i]]$index)]
            names(dates_missing) <- names(stocks[[i]])
            stocks[[i]] <- data.table::copy(dplyr::bind_rows(list(a = dates_missing, b = stocks[[i]])))
            setorder(stocks[[i]], "index")

        }

        n_rows <- unlist(lapply(stocks, nrow))
        keep <- which(n_rows == getmode(n_rows))

        stocks <- stocks[keep]
        SPY_tickers <- names(stocks)			
        SPY_table <- fread(SPY_table_path)[Symbol %in% SPY_tickers]

        #Build Open, High, Low, Close, Volume, Adjusted frames
        stocks[-1] <- lapply(stocks[-1], function(x){if("index" %in% names(x)){x[, index := NULL]}})
        stocks <- data.table::copy(dplyr::bind_cols(stocks))

        varnames <- c("Open", "High", "Low", "Close", "Volume", "Adjusted")
        column_names <- names(stocks)			

        quotes_to_append <- lapply(as.list(varnames), function(x){stocks[, c(1, which(grepl(x, column_names))), with = FALSE]})
        names(quotes_to_append) <- varnames			
        for(i in 1:length(quotes)){

            to_numeric <- names(quotes[[i]])[-1]
            quotes[[i]][, (to_numeric) := lapply(.SD, as.numeric), .SDcols = to_numeric]

            to_numeric <- names(quotes_to_append[[i]])[-1]
            quotes_to_append[[i]][, (to_numeric) := lapply(.SD, as.numeric), .SDcols = to_numeric]

            names(quotes_to_append[[i]]) <- c("Date", SPY_tickers)
            quotes_to_append[[i]][, Date := as.Date(Date)]

            quotes[[i]] <- dplyr::bind_rows(quotes[[i]][!(Date %in% quotes_to_append[[i]]$Date)],
                                            quotes_to_append[[i]])

            filename <- paste("SPY_", stringr::str_to_lower(names(quotes)[i]), ".csv", sep = "")
            path_save <- paste(SPY_quotes_path, filename, sep = "/")
            write.csv(quotes[[i]], path_save, row.names = FALSE)

        }

        n_dates_updated <- nrow(quotes_to_append$Open) - 1
        message(paste("Successfully added", n_dates_updated, "days of historical data for", length(SPY_tickers), "stocks from the SPDR S&P 500 ETF Trust."))


    } else {

        message("Data is already up to date or markets are currently closed. Nothing to query.")

    }

}
