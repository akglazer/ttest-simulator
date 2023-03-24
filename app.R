#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# load libraries
library(shiny)
library(tidyverse)
library(gridExtra)

### functions ###
# simulate random values from a given distribution
simulate_random_values <- function(n = 10, mean = 0, sd = 1, 
                                   rate = 1, shape = 1, pop_vect = c(), distribution = "normal"){
    if(distribution == "normal"){
        x <- rnorm(n, mean, sd)
        pop_mean <- mean
        plot <- ggplot(data.frame(x = c(mean - 3*sd, mean + 3*sd)), aes(x = x)) +
            stat_function(fun = dnorm, args = list(mean = mean, sd = sd)) 
    } else if (distribution == "exponential"){
        x <- rexp(n, rate)
        pop_mean <- 1 / rate
        plot <- ggplot(data.frame(x = c(0, 1/rate + 3*(1/rate^2))), 
                       aes(x = x)) +
            stat_function(fun = dexp, args = list(rate = rate)) 
    } else if (distribution == "gamma") {
        x <- rgamma(n, shape, rate)
        pop_mean <- shape / rate
        plot <- ggplot(data.frame(x = c(0, shape / rate + 3*(shape/rate^2))), 
                       aes(x = x)) +
            stat_function(fun = dgamma, args = list(shape = shape, rate = rate)) 
    } else if (distribution == "custom") {
        x <- sample(pop_vect, n, replace = TRUE)
        pop_mean <- mean(pop_vect)
        plot <- ggplot(data.frame(x = pop_vect), aes(x = x)) +
            geom_histogram()
    }
    
    return(list(sample = x, pop_mean = pop_mean, plot = plot))
}

# simulate null hypothesis rejection rates
simulate_error <- function(n1 = 10, n2 = 10, mean1 = 0, mean2 = 0, sd1 = 1, sd2 = 1, 
                           rate1 = 1, rate2 = 1, shape1 = 1, shape2 = 1, 
                           pop_vect1 = c(), pop_vect2 = c(),
                           dist1 = "normal", dist2 = "normal", reps = 100, alpha = 0.05,
                           add_perm = FALSE, perm_reps = 1000){
    # initialize p-value vector
    p_values <- rep(0, reps)
    # if add permutation test, initialize permutation p-value vector
    if(add_perm){
        t_stats <- rep(0, perm_reps)
        perm_p_values <- rep(0, reps)
    }
    # if custom, run potential R code to set population vector
    if(length(pop_vect1) > 0){
        pop1 <- pop_vect1
    }
    if(length(pop_vect2) > 0){
        pop2 <- pop_vect2
    }
    # if population 2 is same as population 1 set all parameters the same
    if(dist2 == "same_pop1"){
        mean2 <- mean1
        sd2 <- sd1
        rate2 <- rate1
        shape2 <- shape1
        pop2 <- pop1
        dist2 <- dist1
    }
    
    for(i in 1:reps){
        x_values <- simulate_random_values(n = n1, mean = mean1, sd = sd1, rate = rate1, shape = shape1, 
                                    pop_vect = pop1, distribution = dist1)
        x <- x_values$sample
        x_pop_mean <- x_values$pop_mean
        x_plot <- x_values$plot + labs(title = "Distribution of population 1")
        
        y_values <- simulate_random_values(n = n2, mean = mean2, sd = sd2, rate = rate2, shape = shape2, 
                                    pop_vect = pop2, distribution = dist2)
        y <- y_values$sample
        y_pop_mean <- y_values$pop_mean
        y_plot <- y_values$plot + labs(title = "Distribution of population 2")
        
        # calculate p-value
        ttest_values <- t.test(x, y)
        p_values[i] <- ttest_values$p.value
            
        # calculate permutation p-value, if add_perm
        if(add_perm){
            nx <- length(x)
            z <- c(x, y)
            for(j in 1:perm_reps){
                og_value <- ttest_values$statistic
                z <- sample(z, length(z), replace = F)
                t_stats[j] <- t.test(z[1:nx], z[(nx+1):length(z)])$statistic
            }
            perm_p_values[i] <- (sum(t_stats >= og_value) + 1) / (perm_reps + 1)
        }
    }
    
    type_i_error <- mean(p_values < alpha)
    if(add_perm){
        perm_reject_null <- mean(perm_p_values < alpha)   
    } else {
        perm_reject_null <- "No permutation test p-value computed"
    }
    
    return(list(type_i_error = type_i_error, perm_reject_null = perm_reject_null, 
                x_pop_mean = x_pop_mean, y_pop_mean = y_pop_mean, 
                x_plot = x_plot, y_plot = y_plot))
}

### Define UI for application ###
ui <- fluidPage(

    # Application title
    titlePanel("t-test simulator"),
         "This shiny app simulates type I error rate and power for two samples from 
         populations with distributions set below. Type I error rate / power is calculated
over the number of simulations specified below. If you choose a custom population, write population
as a vector as you would in R syntax. Calculating permutation test values can add 
    several minutes to computation time.",
    sidebarLayout(
        sidebarPanel(
            selectInput(
                "dist1", "Distribution 1",
                c(Normal = "normal",
                  Exponential = "exponential",
                  Gamma = "gamma",
                  Custom = "custom")),
            sliderInput("n1",
                        "Sample Size 1:",
                        min = 1,
                        max = 50,
                        value = 30),
            conditionalPanel(
                condition = "input.dist1 == 'normal'",
                numericInput("mean1", "Mean of Population 1:", 0, min = -100, max = 100),
                numericInput("sd1", "SD of Population 1:", 10, min = 1, max = 100)),
            conditionalPanel(
                condition = "input.dist1 == 'exponential'",
                numericInput("rate1", "Rate of Population 1:", 10, min = 1, max = 100)),
            conditionalPanel(
                condition = "input.dist1 == 'gamma'",
                numericInput("shape1", "Shape of Population 1:", 10, min = 1, max = 100),
                numericInput("rate1", "Rate of Population 1:", 10, min = 1, max = 100)),
            conditionalPanel(
                condition = "input.dist1 == 'custom'",
                textInput("customPop1", "Custom Population (Written in R Syntax):", value = "c(rnorm(1000), rnorm(2000, 10))")),
            selectInput(
                "dist2", "Distribution 2",
                c(Normal = "normal",
                  Exponential = "exponential",
                  Gamma = "gamma",
                  Custom = "custom",
                  'Same as Population 1'  = "same_pop1")),
            sliderInput("n2",
                        "Sample Size 2:",
                        min = 1,
                        max = 50,
                        value = 30),
            conditionalPanel(
                condition = "input.dist2 == 'normal'",
                numericInput("mean2", "Mean of Population 2:", 0, min = -100, max = 100),
                numericInput("sd2", "SD of Population 2:", 10, min = 1, max = 100)),
            conditionalPanel(
                condition = "input.dist2 == 'exponential'",
                numericInput("rate2", "Rate of Population 2:", 10, min = 1, max = 100)),
            conditionalPanel(
                condition = "input.dist2 == 'gamma'",
                numericInput("shape2", "Shape of Population 2:", 10, min = 1, max = 100),
                numericInput("rate2", "Rate of Population 2:", 10, min = 1, max = 100)),
            conditionalPanel(
                condition = "input.dist2 == 'custom'",
                textInput("customPop2", "Custom Population (Written in R Syntax):", value = "c(rnorm(1000), rnorm(2000, 10))")),
            numericInput("reps", "Number of simulations:", 100, min = 1, max = 10000),
            numericInput("alpha", "Significance level:", 0.05, min = 0, max = 1),
            
            checkboxInput("add_perm", "Calculate permutation t-test", value = FALSE, width = NULL),
            conditionalPanel(
                condition = "input.add_perm == 1",
                numericInput("perm_reps", "Number of permutations for each permutation test:", 100, min = 10, max = 10000)),
            
            actionButton("simulate", "Run simulation!")
        ),
        

        # show the simulation results
        mainPanel(style="text-align:center; 
                  font-size:20px",
           div(verbatimTextOutput("simResults")),
           plotOutput("distPlot1"),
           plotOutput("distPlot2")
        )
    )
)

### Define server logic ###
server <- function(input, output) {

    n1 <- eventReactive(input$simulate, {input$n1})
    n2 <- eventReactive(input$simulate, {input$n2})
    mean1 <- eventReactive(input$simulate, {input$mean1})
    mean2 <- eventReactive(input$simulate, {input$mean2})
    sd1 <- eventReactive(input$simulate, {input$sd1})
    sd2 <- eventReactive(input$simulate, {input$sd2})
    rate1 <- eventReactive(input$simulate, {input$rate1})
    rate2 <- eventReactive(input$simulate, {input$rate2})
    shape1 <- eventReactive(input$simulate, {input$shape1})
    shape2 <- eventReactive(input$simulate, {input$shape2})
    customPop1 <- eventReactive(input$simulate, {input$customPop1})
    customPop2 <- eventReactive(input$simulate, {input$customPop2})
    dist1 <- eventReactive(input$simulate, {input$dist1})
    dist2 <- eventReactive(input$simulate, {input$dist2})
    add_perm <- eventReactive(input$simulate, {input$add_perm})
    perm_reps <- eventReactive(input$simulate, {input$perm_reps})
    reps <- eventReactive(input$simulate, {input$reps})
    alpha <- eventReactive(input$simulate, {input$alpha})

    output$simResults <- renderText({
        results <- simulate_error(n1 = n1(), n2 = n2(), 
                                  mean1 = mean1(), mean2 = mean2(),
                                  sd1 = sd1(), sd2 = sd2(), 
                                  rate1 = rate1(), rate2 = rate2(),
                                  shape1 = shape1(), shape2 = shape2(),
                                  pop_vect1 = eval(parse(text = customPop1())), 
                                  pop_vect2 = eval(parse(text = customPop2())),
                                  dist1 = dist1(), dist2 = dist2(), 
                                  reps = reps(), alpha = alpha(), add_perm = add_perm(),
                                  perm_reps = perm_reps())
        
        # generate text
        paste("Population 1 has a mean of ", results$x_pop_mean, 
              "and population 2 has a mean of ", results$y_pop_mean, "\n", "\n",
            "Your parametric t-test rejection rate is ", "\n", results$type_i_error, "\n",
              "Your permutation t-test rejection rate is ", "\n", results$perm_reject_null)
    })
    
    output$distPlot1 <- renderPlot(
        {
            results <- simulate_error(n1 = n1(), n2 = n2(), 
                                      mean1 = mean1(), mean2 = mean2(),
                                      sd1 = sd1(), sd2 = sd2(), 
                                      rate1 = rate1(), rate2 = rate2(),
                                      shape1 = shape1(), shape2 = shape2(),
                                      pop_vect1 = eval(parse(text = customPop1())), 
                                      pop_vect2 = eval(parse(text = customPop2())),
                                      dist1 = dist1(), dist2 = dist2(), 
                                      reps = 1, add_perm = F)
            
            results$x_plot
        }
    )

    output$distPlot2 <- renderPlot(
        {
            results <- simulate_error(n1 = n1(), n2 = n2(), 
                                      mean1 = mean1(), mean2 = mean2(),
                                      sd1 = sd1(), sd2 = sd2(), 
                                      rate1 = rate1(), rate2 = rate2(),
                                      shape1 = shape1(), shape2 = shape2(),
                                      pop_vect1 = eval(parse(text = customPop1())), 
                                      pop_vect2 = eval(parse(text = customPop2())),
                                      dist1 = dist1(), dist2 = dist2(), 
                                      reps = 1, add_perm = F)
            
            results$y_plot
        }
    )    
}

# Run the application 
shinyApp(ui = ui, server = server)
