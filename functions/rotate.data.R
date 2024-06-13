rotate.data <- function(x, y, angle){

x.pos <- (x * cos(angle * (pi/180))) + (y * sin((90-angle) * (pi/180)))
y.pos <- (-x * sin(angle * (pi/180))) + (y * cos((90-angle) * (pi/180)))

return(cbind(x, y , x.pos, y.pos))
}