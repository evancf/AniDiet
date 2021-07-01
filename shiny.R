
library("raster")
library("sf")
library("tidyverse")

continent_raster <- raster("Andermann/continent_shapes/raster_sum_all_regions.tif")

continent_proj <- continent_raster %>% 
  as("SpatialPixelsDataFrame") %>% 
  spTransform(CRS('+proj=longlat +datum=WGS84')) %>% 
  st_as_sf() %>% 
  mutate(cell = 1:nrow(.),
         terrestrial = raster_sum_all_regions > 0) %>% 
  dplyr::select(-raster_sum_all_regions)
cell_coords <- continent_proj %>% st_coordinates() %>% 
  as.data.frame() %>% 
  mutate(cell = continent_proj$cell)
colnames(cell_coords) <- c("lng", "lat", "cell")

cell_coords$cell[which(abs(cell_coords$lng - input$MAPID_click$lng) == 
                         min(abs(cell_coords$lng - input$MAPID_click$lng)) & 
                         abs(cell_coords$lat - input$MAPID_click$lat) == 
                         min(abs(cell_coords$lat - input$MAPID_click$lat)))]

# tibble(cell = which(raster::values(continent_raster) != 0),
#        continent = "any")
# 
# continent_raster <- as(continent_raster, "SpatialPixelsDataFrame") %>% 
#   as.data.frame(continent_raster) %>% 
#   mutate(raster_sum_all_regions = ifelse(raster_sum_all_regions > 0, "green", "white"))
# 
# ggplot(continent_raster, aes(x = x, y = y)) + 
#   geom_tile(fill = continent_raster$raster_sum_all_regions)




load("hindcast_webs.RData")


library(shiny)
library(visNetwork)
library(htmlwidgets)
library(leaflet)


# Server -----------------
server <- function(input, output){

  ### Map
  output$mymap <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Esri.WorldImagery,
                       options = providerTileOptions(minZoom = 2, maxZoom = 5)
      )# %>% addMarkers(data = points())
  })
  
  ### Data controls
  
  i <- 9802 #23325 # 120 #
  
  ### Network
  
  output$network_proxy_nodes <- renderVisNetwork({
    nodes <- tibble(id = web_pres_nat[[i]] %>% unlist() %>% unique() %>% sort())
    edges <- data.frame(from = web_pres_nat[[i]][,1], 
                        to = web_pres_nat[[i]][,2])
    
    visNetwork(nodes, edges) %>%
      visExport() %>%
      visEdges(arrows = "to") %>% 
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
                 nodesIdSelection = list(enabled = TRUE)) %>% 
      visPhysics(maxVelocity = 1) %>% 
      visInteraction(zoom = F)# %>% # hover = TRUE, 
      # visEvents(click = "function(nodes) {
      #   Shiny.onInputChange('current_node_id', nodes);
      # ;}")
    
  })
  
  output$view_id <- renderText({
    paste("Current node selection : ", input$network_id_selected)
  })
  
  output$code_network_id <- renderText({
    '
  visNetwork(nodes, edges, main = "Title", submain = "Subtitle") %>%
    visExport() %>%
    visOptions(highlightNearest = TRUE,
      nodesIdSelection = list(enabled = TRUE, selected = "1"))
 '
  })
  
}

# UI -------------------
ui <- fluidPage(
  
  ### Title and instructions
  titlePanel("Food webs without people"),
  #h4("How have extinctions and range changes affected mammal predator-prey interactions?"),
  h5("Click on the map to see how today's food web may have looked without Pleistocene extinctions."),
  
  ##### Top row
  
  ### Map
  fluidRow(
    leafletOutput("mymap"),
  ),
  
  
  ##### Bottom row
  fluidRow(
    
    ### Data controls
    column(
      width = 2,
      selectInput('dataset', 
                  'Choose a dataset:', 
                  choices = c("Natural" = "1", "Current" = "2")),
    ),
    
    ### Network
    column(
      width = 10,
      visNetworkOutput("network_proxy_nodes", height = "800px")
    )
  )#,

  #fluidRow(verbatimTextOutput("code_network_id"))
  
  
)
shinyApp(ui = ui, server = server)
