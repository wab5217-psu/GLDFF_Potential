

typedef struct cell{
  double lat[4];
  double lon[4];
  double center_lat;
  double center_lon;
  int lat_indx;
  int lon_indx;
  int edge;
}CELL;




typedef struct Neighbor{
  int left[2];
  int right[2];
  int lower[3];
  int upper[3];
}NEIGHBOR;


CELL* grid;
NEIGHBOR* neighbors;


/**** In your main ****

  grid=(CELL *)calloc(ngrid,sizeof(CELL));
  neighbors=calloc(ngrid,sizeof(struct Neighbor));
  find_neighbors();


*****************/



int cell_in_out(double lat, double lon, CELL cell )
{
  double dlon;
  double in_lon;
  int j,c;

  in_lon=lon;
  if(in_lon == 360) in_lon=359.8;
  if(in_lon == 0) in_lon=.2;
  
  if( (lat<cell.lat[0]) || (lat>cell.lat[2]) )return(0);
  dlon=cell.lon[1]-cell.lon[0];
    
  if( (in_lon>=cell.lon[0])&&(in_lon<cell.lon[1]) )
    return(1);


  if( ((in_lon-360)>=cell.lon[0])&&((in_lon-360)<cell.lon[1]) )
    return(1);

  if( ((in_lon+360)>=cell.lon[0])&&((in_lon+360)<cell.lon[1]) )
    return(1);
  
  if( cell.lon[1]<cell.lon[0] ){
    if( (in_lon>=cell.lon[0]-360)&&(in_lon<cell.lon[1]) ){
      return(1);
    }
    if( (in_lon>cell.lon[0])&&(in_lon<cell.lon[1]+360) ){
      if( (in_lon>=cell.lon[0]-360)&&(in_lon<cell.lon[1]) ){
	return(1);
      }
    }
  }
  return(0);
}  


void find_neighbors(){
  int jg,ig;
  int jc;
  double lat,lon;
  double lon_r;
  double dlat_l,dlon_l;

  for( jg=0; jg<ngrid; jg++ ){
    /* fprintf(stderr,"find_neighbors %d\n",jg); */
    for( jc=0; jc<3; jc++ )neighbors[jg].lower[jc]=BAD_INT;
    for( jc=0; jc<3; jc++ )neighbors[jg].upper[jc]=BAD_INT;
    for( jc=0; jc<2; jc++ )neighbors[jg].left[jc]=BAD_INT;
    for( jc=0; jc<2; jc++ )neighbors[jg].right[jc]=BAD_INT;
    dlat_l=grid[jg].lat[2]-grid[jg].lat[0];
    dlon_l=fabs(grid[jg].lon[1]-grid[jg].lon[0]);

    if( dlon_l>300 )dlon_l-=360;
    dlon_l=fabs(dlon_l);
    
    /* find lower neighbor of lower-left corner: */
    lat=grid[jg].lat[0]-.1;
    lon=grid[jg].lon[0]+.001;
    for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	neighbors[jg].lower[0]=ig;
	break;
      }
    
    /* step in longitude to find all neighbor cells below current cell */
    ig=neighbors[jg].lower[0];

    for( jc=1; jc<3; jc++ ){
      if( (ig<0) || (ig+jc>ngrid-1) )break;

      lat=grid[ig+jc].lat[3]+.1;
      lon=grid[ig+jc].lon[3]+.001;
      if(cell_in_out(lat,lon,grid[jg])){
	neighbors[jg].lower[jc]=ig+jc;
      }
    }
        
    /* find adjacent right cells */
    lat=grid[jg].lat[0]+.1;
    lon=grid[jg].lon[1]+.1;
    for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	neighbors[jg].right[0]=ig;
	break;
      }
    
    if( dlat_l > (grid[neighbors[jg].right[0]].lat[2]-grid[neighbors[jg].right[0]].lat[0]) ){
      lat=grid[jg].lat[0]+.75*dlat_l;
      lon=grid[jg].lon[1]+.1;
      for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	  neighbors[jg].right[1]=ig;
	  break;
	}
    }    
    
    /* find adjacent left cells */ 
    lat=grid[jg].lat[0]+.1;
    lon=grid[jg].lon[0]-.1;
    for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	neighbors[jg].left[0]=ig;
	break;
      }

    if(neighbors[jg].left[0]<0)fprintf(stderr,"bad neighbor %d %d %6.2f %6.2f\n",jg,neighbors[jg].left[0],grid[jg].center_lat,grid[jg].center_lon);
    /* if( grid[jg].lat[2]!=grid[neighbors[jg].left[0]].lat[2] ){ */
    if( dlat_l > (grid[neighbors[jg].left[0]].lat[2]-grid[neighbors[jg].left[0]].lat[0]) ){
      lat=grid[jg].lat[0]+.75*dlat_l;
      lon=grid[jg].lon[0]-.1;
      for(ig=0; ig<ngrid; ig++)if(cell_in_out(lat,lon,grid[ig])){
	  neighbors[jg].left[1]=ig;
	  break;
	}
    }
    /* fprintf(stderr,"%d lat %5.2f dlat %5.2f  right: %d %d  left: %d %d\n",jg,grid[jg].lat[0],grid[jg].lat[3]-grid[jg].lat[0],neighbors[jg].right[0],neighbors[jg].right[1],neighbors[jg].left[0],neighbors[jg].left[1]); */
  }
}

