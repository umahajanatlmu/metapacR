require(hexSticker)
require(showtext)
font_add_google('Fredoka One', 'myfont')
showtext_auto()
sticker(
  # image
  "man/figures/metabolomeMAP.png", # https://unsplash.com/photos/7wBFsHWQDlk
  s_x=1, # slightly to right to appear centered
  s_y=1,
  s_width=1.5,
  s_height=1.8,
  white_around_sticker = TRUE,

  # package name
  package="metapacR",
  p_size=90,
  p_color = "#e8e8e8", # 00030A 010101 #383838
  p_y = 1,
  p_family = "myfont",
  p_fontface = "bold",

  # Output file
  filename="inst/figures/metapacR_sticker.png",

  # Background colour
  h_fill = "#F0F0F0", # #F0F0F0


  # Border
  # Grey colours: https://www.w3schools.com/colors/colors_shades.asp
  h_color = "#96c5b0",   # 3F4243 7F2B94 3B2691 4238AF
  h_size = 1.5,

  dpi = 1000 # otherwise the final fantasy image quality is not good
);system("open inst/figures/metapacR_sticker.png")
