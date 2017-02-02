%nastaven�� exportu grafu do obrazku - 0 = vypnuto, 1 = zapnuto 
EXPORT = 0;
% 1)
[s, Fs] = audioread ('xkosci00.wav'); s= s';
delka = length(s);
t = delka / Fs;

str = sprintf('Delka: %d vzorku, %d s. Vzorkovaci frekvence je %d Hz',delka,t,Fs);
disp(str)
% 2)
transform = fft(s);
transform = abs(transform(1:8000));
f = (0:delka/2-1)/delka * Fs;

plot(f, transform);xlabel('f [Hz]'); grid;
if EXPORT == 1
print('2', '-dpng');
end

% 3)
indexmax = find(max(transform) == transform);
xmax = f(indexmax);
ymax = transform(indexmax);

str = sprintf('Maximum modulu spektra je na frekvenci %g Hz',xmax); 
disp(str)

% 4)
% zeros
b = [0.2324 -0.4112 0.2324];

%poles 
a = [1 0.2289 0.4662];
if abs(roots(a)) < 1
    disp('Filtr je stabilni');
else
    disp('Filtr je nestabilni');
end
zplane(roots(b),roots(a));
if EXPORT == 1
    print('4','-dpng');
end

% 5)
H = freqz(b,a,256);
f2=(0:255) / 256 * Fs / 2; 
% H= abs(H);
plot(f2,abs(H));
xlabel('f [Hz]'); grid;
if EXPORT == 1
    print('5','-dpng');
end

% 6)
filtered = filter(b,a,s);
filtered = fft(filtered);
filtered = abs(filtered(1:8000));
plot(f, filtered);xlabel('f [Hz]'); grid;
if EXPORT == 1
    print('6','-dpng');
end
    
% 7)
indexmax_f = find(max(filtered) == filtered);
xmax_f = f(indexmax_f);
ymax_f = transform(indexmax_f);

str = sprintf('Maximum modulu spektra filtru je na frekvenci %g Hz',xmax_f); 
disp(str)

% 8)
%obdelnikovy impuls:
obdelnik = [];
sign = 1;
cnt = 0;
for i =1:320
    if cnt == 8
        cnt = 0;
        if sign == 1
            sign = -1;
        else
            sign = 1;
        end
    end
    obdelnik(i) = sign*0.2;
    cnt = cnt +1;
end

obdelnik = fft(obdelnik);
obdelnik = abs(obdelnik(1:160));
f_obd = (0:320/2-1)/320 * Fs;
%stem(s);grid;
%plot(f, transform);xlabel('f [Hz]'); grid; 
%axis([15136-320 15456-320 -inf inf])
% [C1,lag1] = xcorr(s,obdelnik);
% plot(lag1/Fs,abs(C1),'k')
% ylabel('Amplitude')
% grid on

check = [];
ri = 1;
% for i = 0:delka
%     sum = 0;
%     for cnt = 1:length(obdelnik);
%         index = 0;
%         if (cnt+i) >= 1 && (cnt+i) <= delka
%             index = s(cnt+i);
%         end
%         sum = sum + obdelnik(cnt)*index;
%     end
%     check(ri) = abs(sum);
%     ri = ri + 1;
% end

%plot(check);
%plot(s);

% 9)
k = [-50:50];
r =[];
ri = 1;
for i = -50:50
    sum = 0;
    for cnt = 1:delka;
        index = 0;
        if ((cnt+i) >= 1 && (cnt+i) <= delka)
            index = s(cnt+i);
        end
        sum = sum + s(cnt)*index;
    end
    r(ri) = 1/delka* sum;
    ri = ri +1;
end
plot(k,r);grid;
if EXPORT == 1
print('9','-dpng');
end

% 10)
ten_r = find( k == 10);
str = sprintf('Value R[10] is %g', r(ten_r));
disp(str);


% 11)
%interval 1 pro n 
x_axis = [];
for i = 0:40
    x_axis(i+1) = -1 + i*0.05;
end
%interval 2 pro n+10
y_axis = x_axis;

%vysledna array vysledku
res_array = zeros(40);
res_array2 = zeros(40);
%vypocet countu 
for i = 1:40
    for n = 1:delka
        if (s(n) >= x_axis(i) && s(n) < x_axis(i+1))
            for i2 = 1:40
                sn = 0;
                if (n <= delka - 10)
                    sn = s(n+10);
                end
                if (sn >= y_axis(i2) && sn < y_axis(i2+1))
                    res_array(i, abs(i2-41)) = res_array(i, abs(i2-41)) +1;
                end
            end
        end
    end
end

%mam spocitane county, nyni je potreba je podelit 
%delka*40^2 pro hustotu pravdepodobnosti
control_count = 0;
for i = 1:40
    for i2 = 1:40
        control_count = control_count + res_array(i,i2)/(delka);
        res_array2(i,i2) = res_array(i,i2)/((delka)*0.05*0.05);
        res_array(i,i2) = res_array(i,i2)/((delka));
    end
end
x_axis_plot = x_axis;
y_axis_plot = y_axis;
x_axis_plot(41) = [];
y_axis_plot(41) = [];
pcolor(x_axis_plot,y_axis_plot,res_array2);
colormap(jet);
colorbar();
if EXPORT == 1
    print('11','-dpng');
end

% 13)
%indexy nejvetsich pravdepodobnosti
[max_values max_values_i] = max(res_array);

sum2 = 0;
for i=1:40
    middle_value = (x_axis(i) + y_axis(i+1))/2;
    sum2 = middle_value*res_array2(i, max_values_i(i));
end
sum2 = (1/40) * sum2;

str = sprintf('Odhadnuta hodnota R[10] je %g',sum2);
disp(str);
