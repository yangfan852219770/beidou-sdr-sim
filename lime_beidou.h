#ifndef _LIME_BEIDOU_H
#define _LIME_BEIDOU_H

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <lime/LimeSuite.h>

#include "beidou_sim.h"

// 发送频率
#define TX_FREQUENCY 1561098000
#define TX_SAMPLERATE   2046000
#define TX_BANDWIDTH    4092000
// I/Q采样数量
#define NUM_IQ_SAMPLES  (TX_SAMPLERATE / 10)
#define FIFO_LENGTH     (NUM_IQ_SAMPLES * 2)
#define SAMPLES_PER_BUFFER	(32 * 1024)

typedef struct {
	pthread_t thread;
	pthread_mutex_t lock;
	//int error;

	lms_stream_t stream;
	int16_t *buffer;
} tx_t;

typedef struct {
    pthread_t thread;
    pthread_mutex_t lock;
    //int error;

    int ready;
    pthread_cond_t initialization_done;
} beidou_t;

typedef struct {
    tx_t tx;
    beidou_t beidou;

    int status;
    bool finished;
    int16_t *fifo;
    long head;
    long tail;
    size_t sample_length;

    pthread_cond_t fifo_read_ready;
    pthread_cond_t fifo_write_ready;

    double time;
} sim_t;

extern void *beidou_task(void *arg);
extern int is_fifo_write_ready(sim_t *s);
#endif