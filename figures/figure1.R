
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175','#FFEEAD')#颜色设置


DimPlot(bigseu,group.by = "new.ident",cols = my36colors)#cell_type

DimPlot(bigseu,group.by = "patient",cols=c('HCC1'='#E95C59','HCC2'='#E5D2DD','HCC3'='#BD956A',
                                           'HCC4'='#57C3F3','HCC5'='#E95C59','HCC6'='#F1BB72',
                                           'HCC7'='#8C549C','HCC8'='#D6E7A3','HCC9'='#E59CC4'),shuffle = T)