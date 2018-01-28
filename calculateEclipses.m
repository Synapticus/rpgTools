% Eclipse calculator script
% Predicts eclipses and syzygies (3+ body eclipses) between sets of planets 
% Takes a while to plot because it draws every X individually.

% (c) 2018 Sean Patrick 
% per MIT License
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [] = calculateEclipses ()
close all
Epoch = 2400; %How long are we looking for eclipses? In this case, 100 days
% Units are whatever the orbit is in terms of. I.E. if orbital frequency is 1/10, an Epoch of 100 is 10 years.

% List of planets in order of distance from the sun separated by newlines
planetString = 'Abraxas \n Sseras \n Eldrick \n  Nex \n Brachine \n Folon \n Meket \n Zona \n';

% Orbital frequency in 1/(length of a year). Formatted on new lines per orbit ring.
f = [1/1920 1/1920 1/1920 ...
    1/5460 1/5460 ...
    1/1123 1/1123 1/1123];

% Initial phase in radians. Planets with initial phase 0 are in eclipse at t=0.
th = [0 pi/2 3*pi/4 0 pi  pi/4 2*pi/3 0];

% The "shadow" of each planet, that is, how wide is the area it eclipses
% in terms of radian arc behind it? In order of distance from sun.
shadow = [0.3 0.2 0.1 0.2 0.4 0.2 0.3 0.1];

% The number of planets
Planets = length(f);

% Initializing...
eclipseMatrix = zeros(length(f),length(f),100);
eclipseCount = ones(length(f));
s=2;
syzygyMatrix = zeros(10,4);
lastSyzygy = [0 0 0];

% Main orbital loop
for t = 1:Epoch %At all t over the Epoch
    for j = 1:Planets % For each planet
        for k = 1:Planets %See if another planet lines up
       
            %Check if the current phase of two planets line up within their
            %shadow
            
            if (j<k) %Which planet is closer to the sun? ie. Whose shadow matters?
                jkShadow = shadow(j);
            else
                jkShadow = shadow(k);
            end
            
            if (abs(orbitAngle(t,f(j),th(j)) - orbitAngle(t,f(k),th(k))) < jkShadow && ...
               j ~= k && ...
               f(j) ~= f(k)) 
           
                % Store the time of the eclipse and increment the number of
                % eclipses that occurred.
                eclipseMatrix(j,k,eclipseCount(j,k)) = t;
                eclipseCount(j,k) = eclipseCount(j,k) + 1;
                
                
                %Detect syzygies
                for p = 1:Planets %Check if a 3rd planet is in line
                    
                    if (j < p)
                        pjShadow = shadow(j);
                    else
                        pjShadow = shadow(p);
                    end
                    
                    if (k < p)
                        pkShadow = shadow(k);
                    else
                        pkShadow = shadow(p);
                    end
                    
                    if (abs(orbitAngle(t,f(j),th(j)) - orbitAngle(t,f(k),th(k))) < jkShadow && ...
                        abs(orbitAngle(t,f(p),th(p)) - orbitAngle(t,f(k),th(k))) < pkShadow && ...
                        abs(orbitAngle(t,f(j),th(j)) - orbitAngle(t,f(p),th(p))) < pjShadow && ...
                        f(j) ~= f(k) && f(j) ~= f(p) && f(k) ~= f(p) && ...
                        j ~= k && j ~= p && p ~= k && ...
                        ~isempty(setdiff(unique([p,j,k]),lastSyzygy))) %This is to make sure the same event doesn't get recorded twice per t

                        syzygyMatrix(s,:) = [t, p, j, k];
                        lastSyzygy = unique([p,j,k]); %Might need to do an actual search if more orbits are added?
                        s = s+1;

                    end
                end
            
            end %end syzygy check
                       
        end %end k loop
    end %end j loop
    

end %end t loop



eclipseCalendar = zeros(size(eclipseMatrix));

%Trim the zeros out
    for j = 1:length(f)
        for k = 1:length(f)
            eclipseDates = nonzeros(eclipseMatrix(j,k,:));
            if (~isempty(eclipseDates))
                eclipseCalendar(j,k,1:length(nonzeros(eclipseMatrix(j,k,:)))) = nonzeros(eclipseMatrix(j,k,:));
                eclipseCalendar(j,k,length(nonzeros(eclipseMatrix(j,k,:)))+1:end) = zeros(1,length(eclipseMatrix(j,k,:))-length(nonzeros(eclipseMatrix(j,k,:))));
            else
            end
        end
    end
eclipseCalendar(:,:,all(~eclipseCalendar,3)) = [];

%Plotting

%Color gradient used for the set of planets
gradient = 1/Planets:1/Planets:1;
colors = [gradient;zeros(1,length(gradient));1-gradient];

for p1 = 1:8
    figure()
    
    for p2 = 1:8
        
        %Don't waste time plotting within the same orbit. 
        if (f(p1) == f(p2))
            continue
        end
         
                %Trim zeros
                eclipseTimes = squeeze(eclipseCalendar(p1,p2,:));
                
                hold on
                
                %Mark an X at the other planets in an eclipse line
                h =plot(eclipseTimes, p2, 'x');  
                
                %Also mark an X at this planet
                h2 =plot(eclipseTimes, p1, 'x');  
                
                %Set the color of the X
                set(h,'Color',colors(:,p2))
                set(h2,'Color','k')
                
                axis([1 Epoch 0 Planets+1]);
                set(gca, 'XTick', Epoch/10:Epoch/10:Epoch)
                
                ytlblstr = sprintf(planetString);
                ytlbls = regexp(ytlblstr, '\n', 'split');
                ytlbls{end} = 'PLANETS';
                set(gca, 'YTick',1:length(ytlbls), 'YTickLabel',ytlbls)
               
                title(sprintf('Eclipse calendar for %s',ytlbls{p1}));
                xlabel('Hours') %I picked hours for my epoch unit, but it's arbitrary

        % Draw a vertical line where syzygies occurred
        for ind = 1:length(syzygyMatrix(:,1))
            if ( (syzygyMatrix(ind,2) == p1) || ...
                 (syzygyMatrix(ind,3) == p1) || ...
                 (syzygyMatrix(ind,4) == p1) )
                
                line([syzygyMatrix(ind,1), syzygyMatrix(ind,1)],[0 9],'Color',[0.5 0.5 0.5]);   
                
            end
        end
        
    end
end

end

%Find the phase of a planet at time t given year length 1/f and
%initial phase th
function theta = orbitAngle(t,f,th)
theta = mod(2*pi*(f).*t + th,2*pi);
end