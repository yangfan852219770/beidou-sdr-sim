# beidou-sdr-sim


模拟的PRN号为1-5GEO卫星,使用LimeSDR Mini发射。

## 在Linux上安装LimeSDR Mini驱动
sudo add-apt-repository -y ppa:myriadrf/drivers
sudo apt-get update
sudo apt-get install limesuite liblimesuite-dev limesuite-udev limesuite-images
sudo apt-get install soapysdr soapysdr-module-lms7

## 编译运行
make
./test