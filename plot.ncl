load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$GEODIAG_ROOT/geodiag.ncl"
begin
begin
;--- Open file and read in data
  f    = addfile("data_small_rotation.nc","r")
  u  = f->height(0,:,:)                

  wks  = gsn_open_wks("png","plot")    ; send graphics to PNG file
;---Set some resources for color contour lines
  res                  = True
  res@cnMonoLineColor  = False           ; Tells NCL not to draw contour lines in one color
  res@cnLineColors     = span_color_rgba ("NCV_jet",11)   ; NCV_jet has 256 colors; span it to get 11 colors
  res@cnLineThicknessF = 5.0             ; Make lines thicker

  plot = gsn_csm_contour(wks,u,res)

end
end
