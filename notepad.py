#!/usr/bin/env python3
python main.py \
-d /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/ \
-o /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE110530/ \
-m /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE110530/model/ \
--pred /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE110530/m.csv \
-g /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE110530/c.tsv \
-po /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE110530/test_predictions.txt \
-n 100 \
-c 100 \
--steps 100 \
--seed 100 \
--scaling frac \
--prediction_scaling None \
--var_cutoff 0 \
--config GSE110530

42.c100.n1000.s4000.b64.l0.0001

python main.py \
-d /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/ \
-o /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE181034/42.c100.n1000.s4000.b64.l0.0001/ \
-m /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE181034/42.c100.n1000.s4000.b64.l0.0001/model \
--pred /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE181034/m.csv \
-g /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE181034/c.tsv \
-po /mnt/ceph/users/wkoh1/scaden/pipeline/shared_test/Scaden_testset/GSE181034/42.c100.n1000.s4000.b64.l0.0001/test_predictions.txt \
--scaling frac \
--prediction_scaling None \
--var_cutoff 0 \
--seed 10 \
-c 10 \
-n 10 \
--steps 100 \
--batch_size 128 \
--learning_rate 0.0001 \
--config test

