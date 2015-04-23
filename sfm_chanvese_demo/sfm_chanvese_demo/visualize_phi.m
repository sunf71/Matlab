function visualize_phi(phi)
  
  if length(size(phi)) == 2
    [h1 h2] = fat_contour(phi);
    drawnow;
    delete(h1);
    delete(h2);
  end

  if length(size(phi)) == 3
    clf;
    show_vol(phi,'r',1,'lights');
    drawnow;
  end
  
  
  
  
