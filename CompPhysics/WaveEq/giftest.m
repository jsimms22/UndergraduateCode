numfiles = (4000 / 50) + 1;
mydata = cell(1, numfiles);
k = 0;
filename = 'wave.gif';
temp = ones(1,961);
for t = 0:50:4000
      k = k + 1;
      myfilename = sprintf('time%d.dat', t);
      mydata{k} = importdata(myfilename);
        
      dim = mydata{k};
      x = dim(:,1); y = dim(:,2); z = dim(:,3);
      plot3(x,y,z);
      %axis([0 1 0  1 -.5 .5]);
      
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if k == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end