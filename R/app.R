app <- function(...) {
  ui <- ui
  server <- server
  shiny::shinyApp(ui, server, ...)
}

app()
