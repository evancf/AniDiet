# Spatial data -------------------

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


# lat <- 34.5
# lng <- -119
# cell_coords$cell[which(abs(cell_coords$lng - lng) == 
#                          min(abs(cell_coords$lng - lng)) & 
#                          abs(cell_coords$lat - lat) == 
#                          min(abs(cell_coords$lat - lat)))]

# tibble(cell = which(raster::values(continent_raster) != 0),
#        continent = "any")
# 
# continent_raster <- as(continent_raster, "SpatialPixelsDataFrame") %>% 
#   as.data.frame(continent_raster) %>% 
#   mutate(raster_sum_all_regions = ifelse(raster_sum_all_regions > 0, "green", "white"))
# 
# ggplot(continent_raster, aes(x = x, y = y)) + 
#   geom_tile(fill = continent_raster$raster_sum_all_regions)



# Interaction data -------------------------------

load("hindcast_webs.RData")
spp <- unlist(unlist(web_pres_nat)) %>% unique() %>% sort()






# Get wikipedia images -------------

# https://stackoverflow.com/questions/46185030/how-to-download-a-wikipedia-image-for-a-species-page
require(WikipediR); require(rvest)

#titles= vector of page name(s)
#res= desired width in pixels (220 px thumbnail by default)
#savedest= save destination (w terminal '/'); wd by default

getwikipic<-function(titles,res,savedest){
  if(missing(res)){res=220}
  if(missing(savedest)){savedest=NA}
  lapply(titles, function (ttl,...){
    d<-page_info("en","wikipedia",page=ttl,clean_response=T)
    url<-d[[1]]$fullurl
    wikipage<-session(url)
    imginfo<-wikipage %>% html_nodes("tr:nth-child(2) img")
    img.url<- imginfo[1] %>% html_attr("src")
    img.url<-paste0("https:",img.url)
    if(is.na(savedest)){
      savefilename<-paste0(ttl,".jpg")
    }else{savefilename<-paste0(savedest,ttl,".jpg")}
    
    if(res!=220){img.url<-gsub(220,res,img.url)}  
    
    # Skip ones where there isn't a photo available
    if(img.url != "https:"){
      download.file(img.url,savefilename)
      return(paste0("orig.file: ",basename(img.url)))#tell user original filename (or error)
    }
    
  },res,savedest)#End lapply
}

# # Need to get image credits
# picname <- gsub(".*220px-", "", 
#                 img.url)
# 
# paste0("https://en.wikipedia.org/w/api.php?action=query&prop=imageinfo&iiprop=extmetadata&titles=File%3a",
#        picname,
#        "&format=json")

last_downloaded <- which(spp %in% (list.files("./wikipics/") %>% gsub(".jpg", "", .))) %>% tail(1)
if(length(last_downloaded) == 0) last_downloaded <- 1
for(i in last_downloaded:length(spp)){
  # if(i %in% c(175, 216, 278, 335, 337, 339,
  #             366, 378, 380, 490, 518, 526,
  #             529, 584, 640, 655)) next()
  tryCatch(getwikipic(spp[i], 
                      440, 
                      savedest = "./wikipics/")) %>% try()
}
print(i)
spp[i]

# May want to add in some text to these images?
# https://cran.r-project.org/web/packages/magick/vignettes/intro.html




# Write out data for shiny
save(spp, cell_coords, web_current, web_pres_nat, web_no_endangered, file = "global_food_webs_data.RData")




# Shiny app ------------------

library(shiny)
library(visNetwork)
library(htmlwidgets)
library(shinyWidgets) #install.packages("shinyWidgets")
library(leaflet)
library("shinythemes")



# Server -----------------
server <- function(input, output){

  ### Map
  output$mymap <- renderLeaflet({
    leaflet() %>%
      setView(lng = -100, lat = 20, zoom = 2) %>%
      addProviderTiles(providers$Esri.WorldImagery,
                       options = providerTileOptions(minZoom = 1, maxZoom = 4)
      ) %>% addCircleMarkers(lng = -119, 
                      lat = 34.5)
  })
  
  
  #click_layer_id <- NULL
  
  observeEvent(input$mymap_click, {
    
    #if(!is.null(click_layer_id)){
      leafletProxy('mymap') %>% clearMarkers()
    #}
    
    leafletProxy('mymap') %>% addCircleMarkers(lng = input$mymap_click$lng, lat = input$mymap_click$lat)
  })
  
  ### Data controls
  
  #i <- 9802 #23325 # 120 #
  
  clicked_cell <- reactive({
    if(is.null(input$mymap_click$lng)){
      cell <- 11941
    } else{
      cell <- cell_coords$cell[which(abs(cell_coords$lng - input$mymap_click$lng) == 
                                       min(abs(cell_coords$lng - input$mymap_click$lng)) & 
                                       abs(cell_coords$lat - input$mymap_click$lat) == 
                                       min(abs(cell_coords$lat - input$mymap_click$lat)))]
    }
    
    return(cell)
  })
  
  
  # This returns the correct dataset
  datasetInput <- reactive({
    if (input$web_type == "Natural"){
      dataset <- web_pres_nat[[clicked_cell()]]
    }
    else if (input$web_type == "Current"){
      dataset <- web_current[[clicked_cell()]]
    }
    return(dataset)
  })
  
  # https://community.rstudio.com/t/changing-datasets-based-on-select-input-in-r-shiny/67891
  

  
  ### Network
  
  output$network_proxy_nodes <- renderVisNetwork({
    nodes <- tibble(id = datasetInput() %>% unlist() %>% unique() %>% sort())
    edges <- data.frame(from = datasetInput()[,1], 
                        to = datasetInput()[,2])
    if(dim(edges)[2] > 0){
      nodes <- nodes %>% mutate(group = factor(ifelse(nodes$id %in% edges[,1], "yes", "no"), labels = c("no", "yes")))
      nodes <- nodes[order(nodes$group),]
    }
    
    visNetwork(nodes, edges) %>%
      visLayout(randomSeed = 444) %>% 
      #visIgraphLayout() %>% 
      # visEvents(type = "once", startStabilizing = "function() {
      #       this.moveTo({scale:0.4})}") %>%
      visExport() %>%
      visEdges(arrows = "to", smooth = F) %>%  # 
      # visGroups(groupname = "yes", color = list(border = "#2B7CE9", background = "#97C2FC")) %>%
      # visGroups(groupname = "no", color = list(border = "#FFA500", background = "#FFFF00")) %>%
      # #visGroups(groupname = "C", color = list(border = "#FA0A10", background = "#FB7E81")) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
                 nodesIdSelection = list(enabled = TRUE, main = "Species to highlight")) %>% 
      visPhysics(maxVelocity = 2,
                 barnesHut = list(gravitationalConstant = -50000),
                 stabilization = T) %>%  # , stabilization = F
      visInteraction(zoom = T) %>%
      visEvents(click = "function(nodes){
        Shiny.onInputChange('click', nodes.nodes[0]);
        ;}"
      )
    # %>% # hover = TRUE, 
      # visEvents(click = "function(nodes) {
      #   Shiny.onInputChange('current_node_id', nodes);
      # ;}")
    
  })
  
  # output$view_id <- renderText({
  #   paste("Current node selection : ", input$network_id_selected)
  # })
  
 #  output$code_network_id <- renderText({
 #    '
 #  visNetwork(nodes, edges, main = "Title", submain = "Subtitle") %>%
 #    visExport() %>%
 #    visOptions(highlightNearest = TRUE,
 #      nodesIdSelection = list(enabled = TRUE, selected = "1"))
 # '
 #  })
  
  ### Mammal Photos
  
  observe({
    #input$gosel
    if(length(input$click) == 1){
    visNetworkProxy("network_proxy_nodes") %>% visGetSelectedNodes()
    }
  })
  
  # observe({
  #   print(input$network_proxy_nodes_selectedNodes)
  # })
  
  output$mammImage <- renderImage({
    
    if(is.null(input$network_proxy_nodes_selectedNodes)){
      fn <- "./wikipics/Blank.jpg"
    } else{
      fn <- paste0(input$network_proxy_nodes_selectedNodes, ".jpg")

      if(fn %in% list.files("./wikipics/")){
        fn <- paste0("./wikipics/", fn)
      } else{
        fn <- "./wikipics/Unknown.jpg"
      }
    }
    
    list(src = fn, width = "100%")
    
  }, deleteFile = F)
  
  
  
}


# UI -------------------
ui <- fluidPage(theme = shinytheme("lumen"),
                
                ### Title and instructions
                titlePanel("Mammal food webs without Late Pleistocene extinctions"),
                #h4("How have extinctions and range changes affected mammal predator-prey interactions?"),
                h4("Using deep learning to reconstruct how today's food webs would have looked without extinction or habitat loss."),
                
                ##### Top row
                
                ### Map
                # fluidRow(
                #   ,
                # ),
                
                
                ##### Bottom row
                sidebarPanel(collapsable = F,# width = "250px",
                  
                  ### Data controls
                  h5("Show food webs that would occur naturally without mammal extinctions and range changes or that currently occur"),
                  radioGroupButtons("web_type", 
                                    label = NULL,#"Food web type",
                                    c("Natural" = "Natural",
                                      "Current" = "Current")),
                  h5("Click on the map to show food webs from a different location"),
                  leafletOutput("mymap"),
                  br(),
                  h5("Click on a species in the food web to view a photo"),
                  imageOutput("mammImage")
                ),
                  
                  ### Network
                mainPanel(
                  tabsetPanel(
                    tabPanel("World wide (food) web",
                      visNetworkOutput("network_proxy_nodes", 
                                       height = "800px")),
                    tabPanel("FAQ/Methods",
                             h4("What do the links mean?"),
                             p("Links between nodes show predator-prey interactions that are likely to occur among species present in the selected region. Using interaction probabilities modelled using neural networks, we show links where the predicted probability of interaction is over 0.5."),
                             p("The arrows point from consumers to their resources. Only mammals that primarily eat other mammals are shown as predators. Note that consumer interactions can also show scavenging."),
                             h4("Why doesn't a species that I know is present show up in the food web?"),
                             p("Only species where we can estimate that at least one interaction is likely (with interaction probability > 0.5) show up in the food web. So if a species doesn't have any mammal predator (say, a bat), then it won't show up in the food web. Or if we can't be relatively certain (again, probability > 0.5) that a predator consumes any particular prey species, then that predator also won't show up in the food web."), 
                             p("This is similar to how food webs are constructed when observed in the field. Researchers typically report the species that were observed to be involved in predator-prey interactions, but this overlooks species that might be present but were not observed interacting."),
                             h4("[Will add more here!]")),
                    tabPanel("Help",
                             h4("Interacting with the food web"),
                             p("Click on a node to see a photo of the species and highlight the species within two degrees of separation. You can click and drag nodes to rearrange the food web."),
                             p("You can zoom in and out!"),
                             h4("Why doesn't a food web show up for a selected region I selected?"),
                             p("Food webs are only shown for regions where there's at least one predator-prey interaction in the natural scenario. For example, we don't make food web predictions for Hawai'i because there are no native mammals that primarily eat mammals."),
                             h4("[Let me know what could be improved or what is confusing and I can add it here!]"))
                  )
                ),
                
                #fluidRow(verbatimTextOutput("code_network_id"))
                
)
shinyApp(ui = ui, server = server)
