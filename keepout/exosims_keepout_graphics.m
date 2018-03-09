%
% EXOSIMS keepout graphics
% 
% turmon 9/2016
%

close all

% where data lives
% (eventual location)
ROOT_DIR = '/Users/turmon/Mike/Projects/Exosims/src/my-tests/support/keepout/data';
%ROOT_DIR = '/tmp';
ROOT_FN = [ROOT_DIR '/ko-demo'];

fn_ko = [ROOT_FN '-ko.csv'];
fn_tl = [ROOT_FN '-targ.csv'];
fn_tm = [ROOT_FN '-time.csv'];

% tabular representations
ko = readtable(fn_ko);
TL = readtable(fn_tl);
tm = readtable(fn_tm);

% MJD origin in Matlab datenums
t_origin = datenum('17-Nov-1858'); 
t0 = t_origin + tm{1,end};
t1 = t_origin + tm{end,end};
t = t_origin + tm{:,'Time'};

N_targ = size(TL, 1);
N_time = size(tm, 1);

% summary
fprintf('Read %d targets at %d times\n', N_targ, N_time);

%% some plots

% average observability -- mean along times
mu = mean(ko{:,3:end}, 2);

%cmap = jet(1024); cmap = cmap(256:end-256,:); 
cmap = hsv(2048); cmap = brighten(cmap(80:540,:), -0.3); % red...green

%% RA-DEC of all targets
figure
scatter(TL{:,'RA'}, TL{:,'DEC'}, 800, mu, '.');
colormap(cmap)
colorbar
axis([0 360 -90 90])
grid on
title(sprintf('Frequency of Observable Days, %s -- %s', datestr(t0), datestr(t1)), ...
      'FontSize', 16)
xlabel('Right Ascension');
ylabel('Declination');
set(gca, 'FontSize', 12);
export_fig /tmp/ko-freq-plot.pdf
export_fig /tmp/ko-freq-plot.png


%return

%% target-by-target visibility
figure;
[~,inx] = sort(TL{:,'RA'});
imagesc(t, [1:N_targ], 1-ko{inx,3:end});
title('Observable Days: Targets Sorted by Right Ascension', 'FontSize', 16);
datetick('x', 'mmm-yy');
xlabel('Time');
ylabel('Target');
set(gca, 'FontSize', 12);
axis tight;
% not so bright!
colormap(cmap/2); 
% export_fig /tmp/ko-vis-plot.pdf
export_fig /tmp/ko-vis-plot.png



%% sequence
figure;

write_movie = true;
if write_movie,
  % MovieType = 'Uncompressed AVI';
  MovieType = 'Motion JPEG AVI';
  fn_m = '/tmp/ko-animation.avi';
  Movie = VideoWriter(fn_m, MovieType);
  Movie.FrameRate = 12; 
  Movie.Quality = 100; % only for motion jpeg
  open(Movie);
end;

for tinx = 1:size(tm, 1),
  %tm1 = tm{tinx,'time'};
  % visible stars
  vis = (ko{:,tinx+2} > 0);
  clf
  plot(TL{vis,'RA'}, TL{vis,'DEC'}, 'g.', 'MarkerSize', 12)
  hold on
  plot(TL{~vis,'RA'}, TL{~vis,'DEC'}, 'r.', 'MarkerSize', 12)
  tm1 = t_origin + tm{tinx,'Time'};
  title(sprintf('Targets at %s', datestr(tm1)));
  xlabel('Right Ascension');
  ylabel('Declination');
  set(gca, 'FontSize', 12);
  axis([0 360 -90 90])
  grid on;
  set(gcf, 'Position', [1000 963 765 375]);
  drawnow
  if write_movie,
    %writeVideo(Movie, im2frame(frame, cm));
    f = getframe(gcf);
    writeVideo(Movie, f);
  else,
    pause(0.3);
  end;
end

if write_movie,
  close(Movie);
  clear Movie;
  % alas, avi2mp4 is broken
  %opts = struct();
  %status = avi2mp4(fn_m, opts, 'qscale', 4, 'v', '0');
end