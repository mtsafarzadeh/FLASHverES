	    2   k820309              15.0        �c9Z                                                                                                           
       gr_bicgInterface.F90 GR_BICGINTERFACE #         @                                        	                #         @                                        	                #         @                                        	               #BNDTYPES              
                                                        p          p            p                          #         @                                        	               #ISRC    #ISOLN    #POISFACT              
                                                      
                                                      
                                      
      #         @                                   	     	               #IVAR1 
   #IVAR2    #NORM              
                                 
                     
                                                      
                                     
       #         @                                        	               #IVAR    #NORM              
                                                      
                                     
       #         @                                        	               #IRHS    #ILHS    #GCELLFLG              
                                                      
                                                      
                                            #         @                                        	               #ISRC    #ISOLN    #POISFACT    #BCRO    #BCRI    #BCVI    #BCPI    #BCSI    #BCZI    #BCYI    #BC_TYPES    #BC_VALUES                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                   p          p            p                                                                                        
     p          p          p            p          p                          #         @                                   !     	               #IVAR "   #NLAYERS #   #GCELLFLG $             
                                 "                     
                                 #                     
                                 $           #         @                                   %     	               #BCTYPETOAPPLY &   #BCTYPEFROMGRID '   #VARINDEX (   #GRIDDATASTRUCT )   #AXIS *   #FACE +   #IDEST ,                                             &                      
                                 '                     
                                 (                     
                                 )                     
                                 *                     
                                 +                     
                                ,           #         @                                   -     	               #ISOLN .   #ISRC /   #BC_TYPES 0   #BC_VALUES 1             
                                 .                     
                                 /                     
                                 0                       p          p            p                                    
                                 1                   
    p          p          p            p          p                             �   .      fn#fn    �   H       GR_BICGINIT       H       GR_BICGFINALIZE    ^  V       GR_BICGINITSLV (   �  �   a   GR_BICGINITSLV%BNDTYPES    H  k       GR_BICGINITSRC $   �  @   a   GR_BICGINITSRC%ISRC %   �  @   a   GR_BICGINITSRC%ISOLN (   3  @   a   GR_BICGINITSRC%POISFACT    s  h       GR_BICGDOTPROD %   �  @   a   GR_BICGDOTPROD%IVAR1 %     @   a   GR_BICGDOTPROD%IVAR2 $   [  @   a   GR_BICGDOTPROD%NORM    �  \       GR_BICGNORM !   �  @   a   GR_BICGNORM%IVAR !   7  @   a   GR_BICGNORM%NORM    w  j       GR_BICGMULTIAX $   �  @   a   GR_BICGMULTIAX%IRHS $   !  @   a   GR_BICGMULTIAX%ILHS (   a  @   a   GR_BICGMULTIAX%GCELLFLG    �  �       BIPCGSTAB    o  @   a   BIPCGSTAB%ISRC     �  @   a   BIPCGSTAB%ISOLN #   �  @   a   BIPCGSTAB%POISFACT    /  @   a   BIPCGSTAB%BCRO    o  @   a   BIPCGSTAB%BCRI    �  @   a   BIPCGSTAB%BCVI    �  @   a   BIPCGSTAB%BCPI    /	  @   a   BIPCGSTAB%BCSI    o	  @   a   BIPCGSTAB%BCZI    �	  @   a   BIPCGSTAB%BCYI #   �	  �   a   BIPCGSTAB%BC_TYPES $   �
  �   a   BIPCGSTAB%BC_VALUES    7  m       GR_BICGBNDRY "   �  @   a   GR_BICGBNDRY%IVAR %   �  @   a   GR_BICGBNDRY%NLAYERS &   $  @   a   GR_BICGBNDRY%GCELLFLG !   d  �       GR_BICGMAPBCTYPE /     @   a   GR_BICGMAPBCTYPE%BCTYPETOAPPLY 0   T  @   a   GR_BICGMAPBCTYPE%BCTYPEFROMGRID *   �  @   a   GR_BICGMAPBCTYPE%VARINDEX 0   �  @   a   GR_BICGMAPBCTYPE%GRIDDATASTRUCT &     @   a   GR_BICGMAPBCTYPE%AXIS &   T  @   a   GR_BICGMAPBCTYPE%FACE '   �  @   a   GR_BICGMAPBCTYPE%IDEST "   �  z       GRID_SOLVEPRECOND (   N  @   a   GRID_SOLVEPRECOND%ISOLN '   �  @   a   GRID_SOLVEPRECOND%ISRC +   �  �   a   GRID_SOLVEPRECOND%BC_TYPES ,   b  �   a   GRID_SOLVEPRECOND%BC_VALUES 