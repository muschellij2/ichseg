library(hexSticker)
library(desc)
desc = desc::description$new()
package = desc$get("Package")
url = desc$get("BugReports")
url = sub("/issues.*", "", url)
# outline = "#0caa41"
outline = "#009E73"
# background = "#0caa41"
# background = "#7142f4"
background = "black"
# p_color = "black"
sticker("icon.png",
        package = package,
        h_fill = background,
        h_color = outline,
        s_width = 0.4,
        s_height = 0.4,
        s_x = 1,
        url = url,
        u_color = "white",
        u_size = 1.25,
        filename = "sticker.png")


usethis::use_build_ignore(
  c("icon.png", "sticker.R", "sticker.png"))
usethis::use_logo("sticker.png")