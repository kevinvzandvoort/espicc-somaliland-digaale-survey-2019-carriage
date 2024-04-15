#' Create the sampling flowchart
#' - Total numbers are hardcoded

#############################################################################################################################
#' Section 1. Setup flowchart
#############################################################################################################################

labels = list(list(label = "Visited all 894 shelters in Digaale",
                   position = "left", type="in"),
              list(label = "No individuals present in\n405 shelters (on all visits)",
                   position = "right", type="out"),
              list(label = "25 households declined consent",
                   position = "right", type="out"),
              list(label = "Enrolled 464 households in\nhousehold survey",
                   position = "left", type="in"),
              list(label = "Collected demographic data for all\n2,049 individuals living in these\nhouseholds",
                   position = "left", type="in"),
              list(label = "Sampled 531 individuals for\ncontact and individual risk factor\nsurvey",
                   position = "left", type="in"),
              list(label = "22 individuals could not be followed up",
                   position = "right", type="out"),
              list(label = "2 individuals declined consent",
                   position = "right", type="out"),
              list(label = "Enrolled 509 participants (from\n462 households) in contact and\nindividual risk factor survey",
                   position = "left", type="in"),
              list(label = "147 individuals could not be\nfollowed up or declined consent",
                   position = "right", type="out"),
              list(label = "Collected 365 nasopharyngeal\nswabs from contact survey\nparticipants",
                   position = "left", type="in"),
              list(label = "Collected additional 88\nnasopharyngeal swabs from\nother household members",
                   position = "left", type="in"),
              list(label = "Collected 453 nasopharyngeal\nswabs in total",
                   position = "right", type="in", same_y = TRUE))

drawArrow = function(from = labels[[1]], to = labels[[2]], type="L"){
  if(type == "L"){
    grid.lines(x = rep(mean(from$x), 2),
               y = c(from$y[2], mean(to$y)))
    grid.lines(x = c(mean(from$x), to$x[1]),
               y = rep(mean(to$y), 2),
               gp = gpar(fill="#000000"),
               arrow = arrow(ends="last", type="closed",
                             length=unit(line_height*h*0.625, "inches")))
  } else if(type == "rt"){
    grid.lines(x = c(from$x[2], mean(to$x)),
               y = rep(mean(from$y), 2))
    grid.lines(x = rep(mean(to$x), 2),
               y = c(mean(from$y), to$y[1]),
               gp = gpar(fill="#000000"),
               arrow = arrow(ends="last", type="closed",
                             length=unit(line_height*h*0.625, "inches")))
  } else if(type == "|"){
    grid.lines(x = c(mean(from$x), mean(to$x)),
               y = c(from$y[2], to$y[1]),
               gp = gpar(fill="#000000"),
               arrow = arrow(ends="last", type="closed",
                             length=unit(line_height*h*0.625, "inches")))
  } else if(type == "-"){
    grid.lines(x = c(from$x[2], to$x[1]),
               y = c(mean(from$y), mean(to$y)),
               gp = gpar(fill="#000000"),
               arrow = arrow(ends="last", type="closed",
                             length=unit(line_height*h*0.625, "inches")))
  } else if(type == "j"){
    grid.lines(x = rep(mean(from$x), 2),
               y = c(from$y[2], mean(to$y)))
    grid.lines(x = c(mean(from$x), to$x[2]),
               y = rep(mean(to$y), 2),
               gp = gpar(fill="#000000"),
               arrow = arrow(ends="last", type="closed",
                             length=unit(line_height*h*0.625, "inches")))
  } else if(type == "]"){
    add = 2*line_height
    grid.lines(x = c(from$x[2], to$x[2]+add),
               y = rep(mean(from$y), 2))
    grid.lines(x = rep(to$x[2]+add, 2),
               y = c(mean(from$y), mean(to$y)))
    grid.lines(x = c(to$x[2]+add, to$x[2]),
               y = rep(mean(to$y), 2),
               gp = gpar(fill="#000000"),
               arrow = arrow(ends="last", type="closed",
                             length=unit(line_height*h*0.625, "inches")))
  } else if(type == "["){
    add = 1*line_height
    grid.lines(x = c(from$x[1], from$x[1]-add),
               y = rep(mean(from$y), 2))
    grid.lines(x = rep(to$x[1]-add, 2),
               y = c(mean(from$y), mean(to$y)))
    grid.lines(x = c(to$x[1]-add, to$x[1]),
               y = rep(mean(to$y), 2),
               gp = gpar(fill="#000000"),
               arrow = arrow(ends="last", type="closed",
                             length=unit(line_height*h*0.625, "inches")))
  }
}

h = 6 #heigt of image, in inch
w = 6 #width of image, in inch
size = 9 #size of text, in pts

#1pt = 1/72 inch

#############################################################################################################################
#' Generates Figure 1: Flowchart of sampling procedure.
#############################################################################################################################

#x11(height=h, width=w)
pdf(file = sprintf("%s/output/%s/figures/pdf/figure1_sampling_flowchart.pdf", analysis_dir, OUTPUT_DIR), height=h, width=w)
grid.newpage()
line_height = (size/72)/h
box_padding = line_height/2
box_margin = line_height*1
box_width = 0.45
shift=line_height
viewport(x = 0.5, y = 0.5, width = (w-(box_margin*h))/w, height = (h-(box_margin*h))/h,
         just = c(0.5, 0.5), name = sprintf("wrapper")) %>% pushViewport()
for(l in 1:length(labels)){
  box = labels[[l]]
  lines = strsplit(box$label, split = "\n")[[1]]
  nlines=length(lines)
  if(l == 1){
    curr_y = 1#-box_margin #first is at the top
    curr_h = nlines*line_height +box_padding*2    
  } else {
    prev_y = curr_y
    prev_h = curr_h
    curr_h = nlines*line_height +box_padding*2
    if(is.null(box$same_y)) curr_y = prev_y-prev_h-box_margin
    else{
      if(curr_h > prev_h) stop("Method not implemented if box is higher")
      curr_y = curr_y - (prev_h - curr_h)/2
    }
  }
  curr_w = box_width * w
  
  box$x = (ifelse(box$position == "left", 0.25, 0.75) - box_width/2) + c(0, box_width) +shift
  box$y = c(curr_y, curr_y-curr_h)
  labels[[l]] = box
  
  viewport(x = ifelse(box$position == "left", 0.25, 0.75) +shift,
           y = curr_y,
           width = box_width,
           height = curr_h,
           just = c(0.5, 1),
           name = sprintf("box%s", l)) %>%
    pushViewport()
  
  grid.roundrect(x=0, y=1, height=1, width=1, just=c(0, 1),
                 gp = gpar(fill=ifelse(box$type == "in", "#FFFFFF", "#AAAAAA")))
  for(j in 1:nlines){
    grid.text(lines[j], x=(box_padding*h)/curr_w,
              y=(nlines+0.5)/(nlines+1) -(j-1)*(line_height)/(curr_h),
              gp=gpar(fontsize=size), hjust=0, vjust=1)  
  }
  upViewport()
}

drawArrow(labels[[1]], labels[[2]])
drawArrow(labels[[1]], labels[[3]])
drawArrow(labels[[1]], labels[[4]], type="|")
drawArrow(labels[[4]], labels[[5]], type="|")
drawArrow(labels[[5]], labels[[6]], type="|")
drawArrow(labels[[6]], labels[[7]], type="L")
drawArrow(labels[[6]], labels[[8]], type="L")
drawArrow(labels[[6]], labels[[9]], type="|")
drawArrow(labels[[9]], labels[[10]], type="L")
drawArrow(labels[[9]], labels[[11]], type="|")
drawArrow(labels[[5]], labels[[12]], type="[")
drawArrow(labels[[11]], labels[[13]], type="rt")
drawArrow(labels[[12]], labels[[13]], type="-")
dev.off()

png(file = sprintf("%s/output/%s/figures/png/figure1_sampling_flowchart.png", analysis_dir, OUTPUT_DIR),
    height=h, width=w, units = "in", res = 300)
grid.newpage()
line_height = (size/72)/h
box_padding = line_height/2
box_margin = line_height*1
box_width = 0.45
shift=line_height
viewport(x = 0.5, y = 0.5, width = (w-(box_margin*h))/w, height = (h-(box_margin*h))/h,
         just = c(0.5, 0.5), name = sprintf("wrapper")) %>% pushViewport()
for(l in 1:length(labels)){
  box = labels[[l]]
  lines = strsplit(box$label, split = "\n")[[1]]
  nlines=length(lines)
  if(l == 1){
    curr_y = 1#-box_margin #first is at the top
    curr_h = nlines*line_height +box_padding*2    
  } else {
    prev_y = curr_y
    prev_h = curr_h
    curr_h = nlines*line_height +box_padding*2
    if(is.null(box$same_y)) curr_y = prev_y-prev_h-box_margin
    else{
      if(curr_h > prev_h) stop("Method not implemented if box is higher")
      curr_y = curr_y - (prev_h - curr_h)/2
    }
  }
  curr_w = box_width * w
  
  box$x = (ifelse(box$position == "left", 0.25, 0.75) - box_width/2) + c(0, box_width) +shift
  box$y = c(curr_y, curr_y-curr_h)
  labels[[l]] = box
  
  viewport(x = ifelse(box$position == "left", 0.25, 0.75) +shift,
           y = curr_y,
           width = box_width,
           height = curr_h,
           just = c(0.5, 1),
           name = sprintf("box%s", l)) %>%
    pushViewport()
  
  grid.roundrect(x=0, y=1, height=1, width=1, just=c(0, 1),
                 gp = gpar(fill=ifelse(box$type == "in", "#FFFFFF", "#AAAAAA")))
  for(j in 1:nlines){
    grid.text(lines[j], x=(box_padding*h)/curr_w,
              y=(nlines+0.5)/(nlines+1) -(j-1)*(line_height)/(curr_h),
              gp=gpar(fontsize=size), hjust=0, vjust=1)  
  }
  upViewport()
}

drawArrow(labels[[1]], labels[[2]])
drawArrow(labels[[1]], labels[[3]])
drawArrow(labels[[1]], labels[[4]], type="|")
drawArrow(labels[[4]], labels[[5]], type="|")
drawArrow(labels[[5]], labels[[6]], type="|")
drawArrow(labels[[6]], labels[[7]], type="L")
drawArrow(labels[[6]], labels[[8]], type="L")
drawArrow(labels[[6]], labels[[9]], type="|")
drawArrow(labels[[9]], labels[[10]], type="L")
drawArrow(labels[[9]], labels[[11]], type="|")
drawArrow(labels[[5]], labels[[12]], type="[")
drawArrow(labels[[11]], labels[[13]], type="rt")
drawArrow(labels[[12]], labels[[13]], type="-")
dev.off()

setEPS(); postscript(sprintf("%s/output/%s/figures/eps/figure1_sampling_flowchart.eps", analysis_dir, OUTPUT_DIR),
                     width = w, height=h)
grid.newpage()
line_height = (size/72)/h
box_padding = line_height/2
box_margin = line_height*1
box_width = 0.45
shift=line_height
viewport(x = 0.5, y = 0.5, width = (w-(box_margin*h))/w, height = (h-(box_margin*h))/h,
         just = c(0.5, 0.5), name = sprintf("wrapper")) %>% pushViewport()
for(l in 1:length(labels)){
  box = labels[[l]]
  lines = strsplit(box$label, split = "\n")[[1]]
  nlines=length(lines)
  if(l == 1){
    curr_y = 1#-box_margin #first is at the top
    curr_h = nlines*line_height +box_padding*2    
  } else {
    prev_y = curr_y
    prev_h = curr_h
    curr_h = nlines*line_height +box_padding*2
    if(is.null(box$same_y)) curr_y = prev_y-prev_h-box_margin
    else{
      if(curr_h > prev_h) stop("Method not implemented if box is higher")
      curr_y = curr_y - (prev_h - curr_h)/2
    }
  }
  curr_w = box_width * w
  
  box$x = (ifelse(box$position == "left", 0.25, 0.75) - box_width/2) + c(0, box_width) +shift
  box$y = c(curr_y, curr_y-curr_h)
  labels[[l]] = box
  
  viewport(x = ifelse(box$position == "left", 0.25, 0.75) +shift,
           y = curr_y,
           width = box_width,
           height = curr_h,
           just = c(0.5, 1),
           name = sprintf("box%s", l)) %>%
    pushViewport()
  
  grid.roundrect(x=0, y=1, height=1, width=1, just=c(0, 1),
                 gp = gpar(fill=ifelse(box$type == "in", "#FFFFFF", "#AAAAAA")))
  for(j in 1:nlines){
    grid.text(lines[j], x=(box_padding*h)/curr_w,
              y=(nlines+0.5)/(nlines+1) -(j-1)*(line_height)/(curr_h),
              gp=gpar(fontsize=size), hjust=0, vjust=1)  
  }
  upViewport()
}

drawArrow(labels[[1]], labels[[2]])
drawArrow(labels[[1]], labels[[3]])
drawArrow(labels[[1]], labels[[4]], type="|")
drawArrow(labels[[4]], labels[[5]], type="|")
drawArrow(labels[[5]], labels[[6]], type="|")
drawArrow(labels[[6]], labels[[7]], type="L")
drawArrow(labels[[6]], labels[[8]], type="L")
drawArrow(labels[[6]], labels[[9]], type="|")
drawArrow(labels[[9]], labels[[10]], type="L")
drawArrow(labels[[9]], labels[[11]], type="|")
drawArrow(labels[[5]], labels[[12]], type="[")
drawArrow(labels[[11]], labels[[13]], type="rt")
drawArrow(labels[[12]], labels[[13]], type="-")
dev.off()

tiff(file = sprintf("%s/output/%s/figures/tiff/figure1_sampling_flowchart.tiff", analysis_dir, OUTPUT_DIR),
     height=h, width=w, units = "in", res = 300)
grid.newpage()
line_height = (size/72)/h
box_padding = line_height/2
box_margin = line_height*1
box_width = 0.45
shift=line_height
viewport(x = 0.5, y = 0.5, width = (w-(box_margin*h))/w, height = (h-(box_margin*h))/h,
         just = c(0.5, 0.5), name = sprintf("wrapper")) %>% pushViewport()
for(l in 1:length(labels)){
  box = labels[[l]]
  lines = strsplit(box$label, split = "\n")[[1]]
  nlines=length(lines)
  if(l == 1){
    curr_y = 1#-box_margin #first is at the top
    curr_h = nlines*line_height +box_padding*2    
  } else {
    prev_y = curr_y
    prev_h = curr_h
    curr_h = nlines*line_height +box_padding*2
    if(is.null(box$same_y)) curr_y = prev_y-prev_h-box_margin
    else{
      if(curr_h > prev_h) stop("Method not implemented if box is higher")
      curr_y = curr_y - (prev_h - curr_h)/2
    }
  }
  curr_w = box_width * w
  
  box$x = (ifelse(box$position == "left", 0.25, 0.75) - box_width/2) + c(0, box_width) +shift
  box$y = c(curr_y, curr_y-curr_h)
  labels[[l]] = box
  
  viewport(x = ifelse(box$position == "left", 0.25, 0.75) +shift,
           y = curr_y,
           width = box_width,
           height = curr_h,
           just = c(0.5, 1),
           name = sprintf("box%s", l)) %>%
    pushViewport()
  
  grid.roundrect(x=0, y=1, height=1, width=1, just=c(0, 1),
                 gp = gpar(fill=ifelse(box$type == "in", "#FFFFFF", "#AAAAAA")))
  for(j in 1:nlines){
    grid.text(lines[j], x=(box_padding*h)/curr_w,
              y=(nlines+0.5)/(nlines+1) -(j-1)*(line_height)/(curr_h),
              gp=gpar(fontsize=size), hjust=0, vjust=1)  
  }
  upViewport()
}

drawArrow(labels[[1]], labels[[2]])
drawArrow(labels[[1]], labels[[3]])
drawArrow(labels[[1]], labels[[4]], type="|")
drawArrow(labels[[4]], labels[[5]], type="|")
drawArrow(labels[[5]], labels[[6]], type="|")
drawArrow(labels[[6]], labels[[7]], type="L")
drawArrow(labels[[6]], labels[[8]], type="L")
drawArrow(labels[[6]], labels[[9]], type="|")
drawArrow(labels[[9]], labels[[10]], type="L")
drawArrow(labels[[9]], labels[[11]], type="|")
drawArrow(labels[[5]], labels[[12]], type="[")
drawArrow(labels[[11]], labels[[13]], type="rt")
drawArrow(labels[[12]], labels[[13]], type="-")
dev.off()