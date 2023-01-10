function [signal2, win] = time_windowing(signal1, Nleft, Nright, plt)
% TIME_WINDOWING - window signal1 by multiplying cosine functions at the
% beginning and at the end.


Nleft = ceil(Nleft);
Nright = floor(Nright);

if nargin<4
    plt = 0;
end

nsig = size(signal1,1);
nbt = size(signal1,2);

it0 = ones(nsig,1);

% Get it0
for jj=1:nbt
    if all(signal1(:,jj)==0)
        it0(:) = jj;
    else
        break
    end
end

it3 = nbt * ones(nsig,1);

% Get it3
for jj=nbt:-1:1
    if all(signal1(:,jj)==0)
        it3(:) = jj;
    else
        break
    end
end

it1 = it0 + Nleft;
it2 = it3 - Nright;

t = kron(1:nbt, ones(nsig, 1));

win = zeros(nsig, nbt);

for ii=1:nsig
    win(ii,it0(ii):it1(ii)-1) = 0.5 * (1-cos(pi./( t(ii,it1(ii)) - t(ii,it0(ii)) ) .* ...
                                    ( t(ii,it0(ii):it1(ii)-1) - t(ii,it0(ii)) ) ) );
    win(ii,it1(ii):it2(ii)-1) = 1;

    win(ii,it2(ii):it3(ii)-1) = 0.5 * (1+cos(pi/( t(ii,it3(ii)) - t(ii,it2(ii)) ) * ...
                                    ( t(ii,it2(ii):it3(ii)-1) - t(ii,it2(ii)) ) ) );
end

%signal2 = kron(win, ones(size(signal1,1),1)) .* signal1;
signal2 = win .* signal1;
winsum = sum(win, 2);

if plt == 1

figure
plot(1:nbt,signal1./max(abs(signal1),[],2))
hold all
plot(1:nbt,win, '--')
xlim([min(it0)-Nleft, max(it3)+Nright])
%plot(1:nbt,max(signal1,[],2)*win)
legend('normalized original signals', 'window')
xlabel('time, samples')
ylabel('amplitude')
title('time windowing')

figure
plot(1:nbt,signal1)
hold all
plot(1:nbt,signal2,'--')
xlim([min(it0)-Nleft, max(it3)+Nright])
xlabel('time, samples')
ylabel('amplitude')
legend('original signal', 'windowed signal')
title('time windowing')

end

function win = time_window(t, it0, it1, it2, it3)
% TIME_WINDOWING -
%

N = length(t);
win = zeros(1, N);

win(1,it0:it1-1) = 0.5 * (1-cos(pi/( t(it1) - t(it0) ) .* ...
                                ( t(it0:it1-1) - t(it0) ) ) );
win(1,it1:it2-1) = 1;

win(1,it2:it3-1) = 0.5 * (1+cos(pi/( t(it3) - t(it2) ) * ...
                                ( t(it2:it3-1) - t(it2) ) ) );
