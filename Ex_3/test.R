library(ggplot2)
library(dplyr)

raw_data <- "
2300    65  1   4400   56  0
750   156  1   3000   65  0
4300   100  1   4000   17  0
2600   134  1   1500    7  0
6000    16  1   9000   16  0
10500   108  1   5300   22  0
10000   121  1  10000    3  0
17000     4  1  19000    4  0 
5400    39  1  27000    2  0
7000   143  1  28000    3  0
9400    56  1  31000    8  0
32000    26  1  26000    4  0
35000    22  1  21000    3  0
100000     1  1  79000   30  0
100000     1  1 100000    4  0 
52000     5  1 100000   43  0
100000    65  1
"

# First, split the raw data into lines
data_lines <- unlist(strsplit(raw_data, "\n"))

# Remove empty lines
data_lines <- data_lines[data_lines != ""]

# Initialize vectors to store data
WBC <- numeric()
Tempo <- numeric()
AG <- integer()

# Loop through each line and extract data
for (line in data_lines) {
  # Split the line into elements
  elements <- unlist(strsplit(line, "\\s+"))
  elements <- elements[elements != ""]
  
  # Check the number of elements in the line
  if (length(elements) == 6) {
    # First three elements are for AG Positive
    WBC <- c(WBC, as.numeric(elements[1]))
    Tempo <- c(Tempo, as.numeric(elements[2]))
    AG <- c(AG, as.integer(elements[3]))
    
    # Next three elements are for AG Negative
    WBC <- c(WBC, as.numeric(elements[4]))
    Tempo <- c(Tempo, as.numeric(elements[5]))
    AG <- c(AG, as.integer(elements[6]))
  } else if (length(elements) == 3) {
    # Only data for AG Positive
    WBC <- c(WBC, as.numeric(elements[1]))
    Tempo <- c(Tempo, as.numeric(elements[2]))
    AG <- c(AG, as.integer(elements[3]))
  } else {
    warning("Unexpected number of elements in line: ", line)
  }
}

# Create the data frame
data_combined <- data.frame(WBC = WBC, Tempo = Tempo, AG = AG)

# Convert AG to factor with labels
data_combined$AG <- factor(data_combined$AG, levels = c(0, 1),
                           labels = c("Negative", "Positive"))

# View the combined data
print(data_combined)