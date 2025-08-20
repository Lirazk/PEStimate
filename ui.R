# TODO (?):
# Poisson/binomial - doesn't work for exclude without history
# Have the plot and calculator modes synchronized - maybe not
# Cookies to remember things between runs?
# Cache results?
# Minify things?

library(shiny)
library(shinyWidgets)
library(bslib)

slider_and_numeric <- function(id, label, min, max, step, value, helptext = "",
                               placement = "bottom",
                               post = NULL) {
  if(length(step) > 1) {
    div(id = id, 
        withMathJax(),
        splitLayout(cellWidths = c("80%", "20%"),
                    sliderTextInput(
                      inputId = id,
                      label = label,
                      choices = step,
                      grid = F, force_edges = T,
                      post = post),
                    helpPopup(NULL, helptext, placement = placement, c("hover"))))
  }
  else {
    div(id = id, 
    splitLayout(cellWidths = c("80%", "20%"),
                sliderInput(
                  inputId = id,
                  label = label,
                  min = min,
                  max = max,
                  step = step,
                  value = value,
                  post = post
                ), 
                helpPopup(NULL, helptext, placement = placement, c("hover"))))
  }
}


helpPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tooltip(bsicons::bs_icon("question-circle"), content, id="tip", placement = placement)
}

disclamir_and_date_text <- HTML("<p align=\"center\"><b><font color = \"red\">The application is intended for research purposes only and is not intended to guide clinical decision making</font></b><br>",
                                paste("Last update", max(format(max(file.info("ui.R")$mtime,
                                                                    file.info("server.R")$mtime,
                                                                    file.info("EmbryoSelection.R")$mtime), "%d-%m-%Y"))),
                                "</p>")


about_panel <- page_fillable(h1("About", align = "center"),
                             # HTML("<p>This application allows users to predict the expected risk reduction when selecting an IVF embryo for transfer based on polygenic risk scores (PRS) for a single disease.</p>"),
                             p("The app predicts the risk reduction when selecting an IVF embryo for transfer based on polygenic risk scores (PRS) for a single disease."),
                             HTML("<p>Use the <strong>Calculator tab</strong> to present various metrics for risk reduction given any set of model parameters.<br>
                               Use the <strong>Plot tab</strong> to generate graphs of the risk reductions vs some of the parameters.</p>"),
                             h2("Model parameters"),
                             h3("Selection strategy"),
                             p("Users should select one of the following embryo selection strategies."),
                             HTML("<ol>
                                   <li>Lowest risk prioritization. Select the embryo with the lowest PRS among all available embryos.</li>
                                   <li>High-risk exclusion. Exclude embryos with a PRS above a “high-risk” cutoff (to be specified), and then select an embryo at random from among the remaining embryos. In case all embryos are high-risk, select a random embryo.</li></ol>"),
                             h3("Embryo parameters"),
                             p("Users should select one of the following models and specify the associated parameters."),
                             HTML("<ol>
                                  <li>A fixed number of live births. There are n euploid embryos, and each will be born when transferred.</li>
                                  <li>A binomial number of live births: There are n euploid embryos, and each will be born with probability p.</li>
                                  <li>A Poisson number of live births: The number of euploid embryos has a Poisson distribution with mean n, and each will be born with probability p.</li></ol>"),
                             h3("Information on family members"),
                             p("Users can specify the disease status and genetic risk of family members."),
                             HTML("<ol>
                                   <li>Number of parents (of the embryos) with a known disease status (0,1, or 2).</li>
                                   <li>Number of parents affected by the disease.</li>
                                   <li>Number of siblings (of the embryos) with a known disease status.</li>
                                   <li>Number of siblings affected by the disease.</li>
                                   <li>The father's PRS percentile.</li>
                                   <li>The mother's PRS percentile.</li>
                                   </ol>"),
                             h3("Disease parameters"),
                             p("Users should specify the following parameters."),
                             HTML("<ol>
                                   <li>PRS accuracy ($R^2$): The proportion of variance in liability to the disease explained by the PRS. This is a commonly used measure of PRS accuracy. $R^2$ is currently between 5-15% for most common polygenic diseases.</li>
                                   <li>The disease prevalence. The proportion of individuals in the (adult) population affected by the disease.</li>
                                   <li>The heritability of the disease ($h^2$). The proportion of variance in the liability of the disease explained by additive genetic factors. (Only required when specifying the disease status of family members.)</li>
                                   </ol>"),
                             h3("Output"),
                             p("The app reports the following outputs."),
                             HTML("<ol>
                                   <li>The baseline risk (i.e., the risk of a randomly selected embryo).</li>
                                   <li>The risk of the embryo selected based on PRS.</li>
                                   <li>The relative risk reduction.</li>
                                   <li>The absolute risk reduction.</li>
                                   <li>The number of patients who would need to screen their embryos to prevent a single future disease case.</li>
                                   </ol>"),
                             h2("References"),
                             HTML("<ol>
                                   <li> <cite>Lencz, T., Backenroth, D., Granot-Hershkovitz, E., Green, A., Gettler, K., Cho, J. H., Weissbrod, O., Zuk, O., & Carmi, S. (2021). Utility of polygenic embryo screening for disease depends on the selection strategy. eLife. <a href=\"https://doi.org/10.7554/elife.64716\" target=\"_blank\" rel=\"noopener noreferrer\">https://doi.org/10.7554/elife.64716</a> </cite> </li>
                                   <li> <cite>Capalbo, A., de Wert, G., Mertes, H., Klausner, L., Coonen, E., Spinella, F., Van de Velde, H., Viville, S., Sermon, K., Vermeulen, N., Lencz, T. & Carmi, S. (2024) Screening embryos for polygenic disease risk: a review of epidemiological, clinical, and ethical considerations. Hum Reprod Update. <a href=\"https://doi.org/10.1093/humupd/dmae012\" target=\"_blank\" rel=\"noopener noreferrer\">https://doi.org/10.1093/humupd/dmae012</a> </cite> </li>
                                   <li> <cite>Klausner, L., Revital, A., Lencz, T. & Carmi, S. (2025). PEStimate: Predicting offspring disease risk after Polygenic Embryo Screening.</cite> </li>
                                   </ol>"),
                             h2("Contact"),
                             p("Please contact us if you find an error or have any suggestion."),
                             HTML("<p>Shai Carmi, <a href=\"mailto: shai.carmi@huji.ac.il\">shai.carmi@huji.ac.il</a></p>"),
                             HTML("<p>Liraz Klausner, <a href=\"mailto: liraz.klausner@mail.huji.ac.il\">liraz.klausner@mail.huji.ac.il</a></p>"),
                             p("Braun School of Public Health, The Hebrew University of Jerusalem"),
                             disclamir_and_date_text)

plot_panel <- page_fillable(
  layout_columns(col_widths = c(4, 8), fill = F,
                 div(
  card(
    conditionalPanel("input.lowestexclude == \"Lowest\"",
                     selectInput("det_random", label = "Is the number of embryos fixed?",
                                 choices = c("Fixed live births", 
                                             "Binomial",
                                             "Poisson"),
                                 width = "50%"),
                     conditionalPanel("input.det_random != \"Fixed live births\"", 
                                      slider_and_numeric("p_lb", "Probability of live birth:",
                                                         0.01, 0.99, 0.01, value=0.3, 
                                                         helptext = "Live birth rate per euploid embryo transfer"))),
    radioButtons("x_var", "Variable for x axis", choiceNames = c("R-squared", "Disease prevalence", "Number of embryos"), 
                 choiceValues = c("r2", "Disease prevalence", "Number of embryos"),
                 selected = "r2",
                 inline = T),
    radioButtons(inputId = "lowestexclude",
                 label = "Choose lowest risk embryo or exclude high risk embroys",
                 choices = c("Lowest", "Exclude"), inline = T
    )),
  card(
    conditionalPanel("input.x_var != \"r2\"",
                     slider_and_numeric("r", "PRS accuracy ($R^2$):", 0.01, 1, NULL, 0.05, 
                                        "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.15.")),
    conditionalPanel("input.x_var != \"Disease prevalence\"",
                     slider_and_numeric("K", "Disease prevalence:", 0.001*100, 0.3*100, 
                                        100*sort(unique(c(seq(0.001, 0.3, 0.001), 
                                                          round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)))), 100*0.001, 
                                        "Fraction of the population with the disease",
                                        post = "%")),
    conditionalPanel("input.x_var != \"Number of embryos\"",
                     slider_and_numeric("N", "Number of embryos:", 2, 10, 1, 5, "The number of embryos/births avaliable for selection.")),
    conditionalPanel(
      condition = "input.lowestexclude == 'Exclude' & input.x_var != 'Percentile'",
      slider_and_numeric("q", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.1, all embryos with PRS at the top 10% of the distribution of PRS in the population will be excluded. If no embryo is avaliable, select one at random.")))
  )),
  card(plotOutput(outputId = "distPlot"))),
  disclamir_and_date_text
)

calc_panel <- page_fillable(
  tags$head(tags$script(("
    Number.prototype.mapLog = function (min, max) {
      const mapped = (this - min) * (Math.log(max) - Math.log(min)) / (max - min) + Math.log(min);
      return Math.exp(mapped);
    }
    $(function() {
    setTimeout(function(){
    
    # $('#K2').data('ionRangeSlider').update({
    #        'prettify': function (n) {
    #        return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
    #        }})

    $('#K_1').data('ionRangeSlider').update({
           'prettify': function (n) {
           return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
           }})
    $('#K_2').data('ionRangeSlider').update({
           'prettify': function (n) {
           return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
           }})
                           }, 2)})"))),
  layout_columns(col_widths = c(4, 8), fill=F, 
                 # layout_column_wrap(heights_equal = "row", width=1,
                 div(
                   card(card_header("Embryo parameters"),
                        conditionalPanel("input.lowestexclude2 == \"Lowest\"",
                                         selectInput("det_random2", label = "Is the number of embryos fixed?",
                                                     choices = c("Fixed live births", 
                                                                 "Binomial",
                                                                 "Poisson"),
                                                     width = "50%"),
                                         conditionalPanel("input.det_random2 != \"Fixed live births\"", 
                                                          slider_and_numeric("p_lb2", "Probability of live birth:",
                                                                             0.01, 0.99, 0.01, value=0.3,
                                                                             helptext = "Live birth rate per euploid embryo transfer"))),
                        slider_and_numeric("N2", "Number of live births:", 2, 10, 1, 5, "The number of embryos/births avaliable for selection."),
                        radioButtons(
                          inputId = "lowestexclude2",
                          label = "Choose lowest risk embryo or exclude high risk embroys",
                          choices = c("Lowest", "Exclude"), inline = T)),
                   card(card_header("Disease parameters"),
                        # slider_and_numeric("K2", "Disease prevalence:", 0.01, 1, 0.01, 0.5, NULL),
                        slider_and_numeric("K2", "Disease prevalence:", 100*0.001, 100*0.3, 
                                           100*sort(unique(c(seq(0.001, 0.3, 0.001), 
                                                             round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)))), 100*0.001, "Fraction of the population with the disease",
                                           post = "%"),
                        slider_and_numeric("r2", "PRS accuracy ($R^2$):", 0.01, 0.99, NULL, 0.05, "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.15. Must be smaller than $h^2$ when conditioning on family disease status."),
                        slider_and_numeric("q2", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.1, all embryos with PRS at the top 10% of the distribution of the PRS in the population will be excluded. If no embryo is avaliable, select one at random.")),
                        slider_and_numeric("h2", "$h^2:$", 0.01, 0.99, 0.01, 0.4, "The heritability of the disease. Only relevant when conditioning on the family disease status."))),
                 card(card_header("Family information", style = "text-align: center;"),
                      layout_column_wrap(width=1/2,
                                         card(card_header("Family disease status",
                                                          helpPopup(NULL, "Do not condition on the disease status of family members who are too young to develop the disease.", "bottom", "hover")),
                                           # conditionalPanel(
                                           # condition = "input.type2 == 'Conditional'",
                                           # ),
                                           # conditionalPanel(
                                           # condition = "input.type2 == 'Family History'",
                                           sliderInput(inputId ="parents_n",
                                                       label = "Number of parents with known disease status:",
                                                       min=0,
                                                       max=2,
                                                       value=0,
                                                       step=1),
                                           sliderInput(inputId = "sick_parents",
                                                       label = "Number of affected parents:",
                                                       min = 0,
                                                       max = 2,
                                                       value = 0,
                                                       step = 1),
                                           sliderInput("siblings_n",
                                                       label = "Number of sibling with known disease status:",
                                                       min = 0,
                                                       max = 20,
                                                       value = 0,
                                                       step = 1),
                                           sliderInput("sick_siblings",
                                                       label = "Number of affected siblings:",
                                                       min = 0,
                                                       max = 20,
                                                       value = 0,
                                                       step = 1)),
                                         div(
                                             card(card_header("Parental polygenic risk score"),
                                               radioButtons(
                                                 inputId = "type2",
                                                 label = "Condition on the parentel polygenic risk score?",
                                                 # choices = c("Risk reduction", "Conditional", "Family History"), inline = T
                                                 # choiceValues = c("Risk reduction", "Conditional", "Family History"),
                                                 # choiceNames = c("No conditioning", "Conditional on the parents' polygenic risk score", "Conditional on family disease status"),
                                                 choiceValues = c("Risk reduction", "Conditional"),
                                                 choiceNames = c("No", "Yes"), inline = T,
                                               ),
                                               slider_and_numeric("qf2", "Father's polygenic risk score percentile:", 0.01, 0.99, 0.01, 0.5, paste("For example, if this parameter equal 0.05, the PRS of the father is at the top 5% of the distribution of the PRS in the population.")),
                                               slider_and_numeric("qm2", "Mother's polygenic risk score percentile:", 0.01, 0.99, 0.01, 0.5, paste("For example, if this parameter equal 0.05, the PRS of the mother is at the top 5% of the distribution of the PRS in the population."))), 
                                             card(card_header("Results"),
                                                  class = "alert alert-info", htmlOutput("summary")))))
  ),
  disclamir_and_date_text
)

calc_two_traits <- div(class = "well", fluidRow(column(4,
                                                       slider_and_numeric("N_2", "Number of embryos:", 2, 10, 1, 5, "The number of embryos available for selection."),
                                                       slider_and_numeric("rho", '$\\rho$, the genetic correlation between the diseases:', -0.99, 0.99, 0.01, 0, "The genetic correlation between the two diseases."),
                                                       slider_and_numeric("samples_2", "Number of monte carlo draws:", 100000, 500000, 1000, 100000, "The number of simulations. Higher number will give a more accurate estimate, but might take longer to run.")),
                                                column(4, 
                                                       slider_and_numeric("r2_1", "PRS accuracy ($R^2 ~ \\text{disease 1}$):", 0.01, 1, 0.001, 0.05, "The proportion of the variance in the liability of the first disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
                                                       slider_and_numeric("r2_2", "PRS accuracy ($R^2 ~ \\text{disease 2}$):", 0.01, 1, 0.001, 0.05, "The proportion of the variance in the liability of the second disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
                                                       fluidRow(column(8, offset = 2, htmlOutput("two_traits"), align = "center"))),
                                                column(4, 
                                                       slider_and_numeric("K_1", "Prevalence of disease 1:", 0.001, 0.3, unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)), 0.001, "How prevalent is the first disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."),
                                                       slider_and_numeric("K_2", "Prevalence of disease 2:", 0.001, 0.3, unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)), 0.001, "How prevalent is the second disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."))),
                       disclamir_and_date_text)

addResourcePath(prefix = "Images", directoryPath = "Images/")

ui <- page_navbar(
  header=withMathJax(tags$head(HTML("<script type='text/x-mathjax-config' >
    MathJax.Hub.Config({
    tex2jax: {inlineMath: [['$','$']]}
    });
    </script>")), shinyjs::useShinyjs()),
  title=tags$a("", uiOutput("logo_ui")),
  nav_panel("About", about_panel),
  nav_panel("Plot", plot_panel),
  nav_panel("Calculator", calc_panel),
  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode")),
)
