#include "lime_beidou.h"

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <lime/LimeSuite.h>

void init_sim(sim_t *s)
{
    pthread_mutex_init(&(s->tx.lock), NULL);
    //s->tx.error = 0;

    pthread_mutex_init(&(s->beidou.lock), NULL);

    s->beidou.ready = 0;
    pthread_cond_init(&(s->beidou.initialization_done), NULL);

    s->status = 0;
    s->head = 0;
    s->tail = 0;
    s->sample_length = 0;

    pthread_cond_init(&(s->fifo_write_ready), NULL);
    pthread_cond_init(&(s->fifo_read_ready), NULL);

    s->time = 0.0;
}

size_t get_sample_length(sim_t *s)
{
    long length;

    length = s->head - s->tail;
    if (length < 0)
        length += FIFO_LENGTH;

    return((size_t)length);
}

size_t fifo_read(int16_t *buffer, size_t samples, sim_t *s)
{
    size_t length;
    size_t samples_remaining;
    int16_t *buffer_current = buffer;

    length = get_sample_length(s);

    if (length < samples)
        samples = length;

    length = samples; // return value

    samples_remaining = FIFO_LENGTH - s->tail;

    if (samples > samples_remaining) {
        memcpy(buffer_current, &(s->fifo[s->tail * 2]), samples_remaining * sizeof(int16_t) * 2);
        s->tail = 0;
        buffer_current += samples_remaining * 2;
        samples -= samples_remaining;
    }

    memcpy(buffer_current, &(s->fifo[s->tail * 2]), samples * sizeof(int16_t) * 2);
    s->tail += (long)samples;
    if (s->tail >= FIFO_LENGTH)
        s->tail -= FIFO_LENGTH;

    return(length);
}

bool is_finished_generation(sim_t *s)
{
    return s->finished;
}

int is_fifo_write_ready(sim_t *s)
{
    int status = 0;

    s->sample_length = get_sample_length(s);
    if (s->sample_length < NUM_IQ_SAMPLES)
        status = 1;

    return(status);
}

void *tx_task(void *arg)
{
    sim_t *s = (sim_t *)arg;
    size_t samples_populated;

    while (1) {
        int16_t *tx_buffer_current = s->tx.buffer;
        unsigned int buffer_samples_remaining = SAMPLES_PER_BUFFER;

        while (buffer_samples_remaining > 0) {

            pthread_mutex_lock(&(s->beidou.lock));
            while (get_sample_length(s) == 0)
            {
                pthread_cond_wait(&(s->fifo_read_ready), &(s->beidou.lock));
            }
//			assert(get_sample_length(s) > 0);

            samples_populated = fifo_read(tx_buffer_current,
                                          buffer_samples_remaining,
                                          s);
            pthread_mutex_unlock(&(s->beidou.lock));

            pthread_cond_signal(&(s->fifo_write_ready));
#if 0
            if (is_fifo_write_ready(s)) {
				/*
				printf("\rTime = %4.1f", s->time);
				s->time += 0.1;
				fflush(stdout);
				*/
			}
			else if (is_finished_generation(s))
			{
				goto out;
			}
#endif
            // Advance the buffer pointer.
            buffer_samples_remaining -= (unsigned int)samples_populated;
            tx_buffer_current += (2 * samples_populated);
        }

        // If there were no errors, transmit the data buffer.
        LMS_SendStream(&s->tx.stream, s->tx.buffer, SAMPLES_PER_BUFFER, NULL, 1000);
        if (is_fifo_write_ready(s)) {
            /*
            printf("\rTime = %4.1f", s->time);
            s->time += 0.1;
            fflush(stdout);
            */
        }
        else if (is_finished_generation(s))
        {
            goto out;
        }

    }
    out:
    return NULL;
}

int start_tx_task(sim_t *s)
{
    int status;

    status = pthread_create(&(s->tx.thread), NULL, tx_task, s);

    return(status);
}

int start_beidou_task(sim_t *s)
{
    int status;

    status = pthread_create(&(s->beidou.thread), NULL, beidou_task, s);

    return(status);
}

int main(){
    sim_t s;
    // 查找设备
	int device_count = LMS_GetDeviceList(NULL);
    // 增益
    double gain = 0.1;

    if (device_count < 1)
	{
		printf("ERROR: No device was found.\n");
		exit(1);
	}
	else if (device_count > 1)
	{
		printf("ERROR: Found more than one device.\n");
		exit(1);
	}
    // 查找设备
    lms_info_str_t *device_list = malloc(sizeof(lms_info_str_t) * device_count);
	device_count = LMS_GetDeviceList(device_list);

    // Initialize simulator
    init_sim(&s);

	// 分配I/Q缓存
	s.tx.buffer = (int16_t *)malloc(SAMPLES_PER_BUFFER * sizeof(int16_t) * 2);
    if (s.tx.buffer == NULL)
    {
        printf("ERROR: 分配tx buffer失败.\n");
        goto out;
    }
    // 分配fifo缓存
    s.fifo = (int16_t *)malloc(FIFO_LENGTH * sizeof(int16_t) * 2);
    if (s.fifo == NULL)
    {
        printf("ERROR: Failed to allocate I/Q sample buffer.\n");
        goto out;
    }

	lms_device_t *device = NULL;
    
    printf("Opening and initializing device...\n");
    // 打开设备
    if(LMS_Open(&device, device_list[0], NULL)){
        printf("ERROR: Failed to open device: %s\n", device_list[0]);
        goto out;
    }
    // 获取设备信息
    const lms_dev_info_t *devinfo =  LMS_GetDeviceInfo(device);
	if (devinfo == NULL)
	{
		printf("ERROR: Failed to read device info: %s\n", LMS_GetLastErrorMessage());
		goto out;
	}
    // 打印设备信息
    printf("deviceName: %s\n", devinfo->deviceName);
    printf("expansionName: %s\n", devinfo->expansionName);
    printf("firmwareVersion: %s\n", devinfo->firmwareVersion);
    printf("hardwareVersion: %s\n", devinfo->hardwareVersion);
    printf("protocolVersion: %s\n", devinfo->protocolVersion);
    printf("gatewareVersion: %s\n", devinfo->gatewareVersion);
    printf("gatewareTargetBoard: %s\n", devinfo->gatewareTargetBoard);

    int limeOversample = 1;
    // 重置设备
    int lmsReset = LMS_Reset(device);
	if (lmsReset)
	{
		printf("ERROR: Failed to reset device: %s\n", LMS_GetLastErrorMessage());
		goto out;
	}

    // 初始化设备
    int lmsInit = LMS_Init(device);
	if (lmsInit)
	{
		printf("ERROR: Failed to linitialize device: %s\n", LMS_GetLastErrorMessage());
		goto out;
	}
    // 选择天线频道
    int32_t channel = 0;
    // 可用天线数量
    int antenna_count = LMS_GetAntennaList(device, LMS_CH_TX, channel, NULL);
	lms_name_t *antenna_name = malloc(sizeof(lms_name_t) * antenna_count);
    // 获取发送端口的带宽
    if (antenna_count > 0)
	{
		int i = 0;
		lms_range_t *antenna_bw = malloc(sizeof(lms_range_t) * antenna_count);
		LMS_GetAntennaList(device, LMS_CH_TX, channel, antenna_name);
		for (i = 0; i < antenna_count; i++)
		{
			LMS_GetAntennaBW(device, LMS_CH_TX, channel, i, antenna_bw + i);
			//printf("Channel %d, antenna [%s] has BW [%lf .. %lf] (step %lf)" "\n", channel, antenna_name[i], antenna_bw[i].min, antenna_bw[i].max, antenna_bw[i].step);
		}
	}
    // 设置增益
    LMS_SetNormalizedGain(device, LMS_CH_TX, channel, gain);
	// 禁用其它频道
	LMS_EnableChannel(device, LMS_CH_TX, 1 - channel, false);
	LMS_EnableChannel(device, LMS_CH_RX, 0, false);
	LMS_EnableChannel(device, LMS_CH_RX, 1, false);
	// 开启要使用的Tx0频道
	LMS_EnableChannel(device, LMS_CH_TX, channel, true);
    int setLOFrequency = LMS_SetLOFrequency(device, LMS_CH_TX, channel, (double)TX_FREQUENCY);
    if (setLOFrequency)
    {
		printf("ERROR: Failed to set TX frequency: %s\n", LMS_GetLastErrorMessage());
		goto out;
	}

    // 设置采样率
    lms_range_t sampleRateRange;
    // 采样率范围
	int getSampleRateRange = LMS_GetSampleRateRange(device, LMS_CH_TX, &sampleRateRange);
    if (getSampleRateRange)
		printf("Warning: Failed to get sample rate range: %s\n", LMS_GetLastErrorMessage());
    int setSampleRate = LMS_SetSampleRate(device, (double)TX_SAMPLERATE, limeOversample); 
	if (setSampleRate)
	{
		printf("ERROR: Failed to set sample rate: %s\n", LMS_GetLastErrorMessage());
		goto out;
	}
    double actualHostSampleRate = 0.0;
	double actualRFSampleRate = 0.0;
    // 获得采样速率
	int getSampleRate = LMS_GetSampleRate(device, LMS_CH_TX, channel, &actualHostSampleRate, &actualRFSampleRate);
	if (getSampleRate)
		printf("Warnig: Failed to get sample rate: %s\n", LMS_GetLastErrorMessage());
	else
		printf("Sample rate: %.1lf Hz (Host) / %.1lf Hz (RF)" "\n", actualHostSampleRate, actualRFSampleRate);
    // 自动校准
	printf("Calibrating...\n");
	int calibrate = LMS_Calibrate(device, LMS_CH_TX, channel, (double)TX_BANDWIDTH, 0);
	if (calibrate)
		printf("Warning: Failed to calibrate device: %s\n", LMS_GetLastErrorMessage());


    // 设置Tx发送流
	printf("Setup TX stream...\n");
	s.tx.stream.channel = channel;
	s.tx.stream.fifoSize = 1024 * 1024;
    s.tx.stream.throughputVsLatency = 0.5;
	s.tx.stream.isTx = true;
	s.tx.stream.dataFmt = LMS_FMT_I12;
	int setupStream = LMS_SetupStream(device, &s.tx.stream);
    if (setupStream)
	{
		printf("ERROR: Failed to setup TX stream: %s\n", LMS_GetLastErrorMessage());
		goto out;
	}
    // 启动Tx流
    LMS_StartStream(&s.tx.stream);
    // 发送北斗数据

    // Start beidou task.
    s.status = start_beidou_task(&s);
    if (s.status < 0) {
        fprintf(stderr, "Failed to start BEIDOU task.\n");
        goto out;
    }
    else
        printf("Creating BEIDOU task...\n");

    // Wait until biedou task is initialized
    pthread_mutex_lock(&(s.tx.lock));
    while (!s.beidou.ready)
        pthread_cond_wait(&(s.beidou.initialization_done), &(s.tx.lock));
    pthread_mutex_unlock(&(s.tx.lock));

    // Fillfull the FIFO.
    if (is_fifo_write_ready(&s))
        pthread_cond_signal(&(s.fifo_write_ready));

    // Start TX task
    s.status = start_tx_task(&s);
    if (s.status < 0) {
        fprintf(stderr, "Failed to start TX task.\n");
        goto out;
    }
    else
        printf("Creating TX task...\n");

    // Wainting for TX task to complete.
    pthread_join(s.tx.thread, NULL);
    printf("\nDone!\n");

out:
    // 关闭Tx流
    LMS_StopStream(&s.tx.stream);
    LMS_DestroyStream(device, &s.tx.stream);
    // 关闭设备
    printf("Closing device...\n");
	LMS_EnableChannel(device, LMS_CH_TX, channel, false);
	LMS_Close(device);

    return 0;
}