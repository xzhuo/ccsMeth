# -*- coding: utf-8 -*-
import argparse
import os
import re
import sys
import time

import numpy as np
import torch
import torch.nn as nn
from sklearn import metrics
from torch.optim.lr_scheduler import StepLR

from .dataloader import FeaData
from .dataloader import FeaData2
from .dataloader import FeaData2s
from .dataloader import clear_linecache

from .models import ModelRNN
from .models import ModelAttRNN
from .models import ModelResNet18
from .models import ModelTransEncoder
from .models import ModelAttRNN2s

from .utils.constants_torch import use_cuda
from .utils.process_utils import display_args
from .utils.process_utils import str2bool


def train(args):
    total_start = time.time()
    torch.manual_seed(args.tseed)
    torch.cuda.manual_seed(args.tseed)

    print("[train]start..")
    if use_cuda:
        print("GPU is available!")
    else:
        print("GPU is not available!")

    print("reading data..")
    if args.model_type in {"bilstm", "bigru", "attbilstm", "attbigru",
                           "transencoder", }:
        train_dataset = FeaData(args.train_file)
        train_loader = torch.utils.data.DataLoader(dataset=train_dataset,
                                                   batch_size=args.batch_size,
                                                   shuffle=True)

        valid_dataset = FeaData(args.valid_file)
        valid_loader = torch.utils.data.DataLoader(dataset=valid_dataset,
                                                   batch_size=args.batch_size,
                                                   shuffle=False)
    elif args.model_type in {"resnet18", }:
        train_dataset = FeaData2(args.train_file)
        train_loader = torch.utils.data.DataLoader(dataset=train_dataset,
                                                   batch_size=args.batch_size,
                                                   shuffle=True)

        valid_dataset = FeaData2(args.valid_file)
        valid_loader = torch.utils.data.DataLoader(dataset=valid_dataset,
                                                   batch_size=args.batch_size,
                                                   shuffle=False)
    elif args.model_type in {"attbigru2s", }:
        train_dataset = FeaData2s(args.train_file)
        train_loader = torch.utils.data.DataLoader(dataset=train_dataset,
                                                   batch_size=args.batch_size,
                                                   shuffle=True)

        valid_dataset = FeaData2s(args.valid_file)
        valid_loader = torch.utils.data.DataLoader(dataset=valid_dataset,
                                                   batch_size=args.batch_size,
                                                   shuffle=False)
    else:
        raise ValueError("model_type not right!")

    model_dir = args.model_dir
    model_regex = re.compile(r"" + args.model_type + "\.b\d+_epoch\d+\.ckpt*")
    if model_dir != "/":
        model_dir = os.path.abspath(model_dir).rstrip("/")
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)
        else:
            for mfile in os.listdir(model_dir):
                if model_regex.match(mfile):
                    os.remove(model_dir + "/" + mfile)
        model_dir += "/"

    if args.model_type in {"bilstm", "bigru", }:
        model = ModelRNN(args.seq_len, args.layer_rnn, args.class_num,
                         args.dropout_rate, args.hid_rnn,
                         args.n_vocab, args.n_embed,
                         is_stds=str2bool(args.is_stds),
                         model_type=args.model_type)
    elif args.model_type in {"attbilstm", "attbigru", }:
        model = ModelAttRNN(args.seq_len, args.layer_rnn, args.class_num,
                            args.dropout_rate, args.hid_rnn,
                            args.n_vocab, args.n_embed,
                            is_stds=str2bool(args.is_stds),
                            model_type=args.model_type)
    elif args.model_type in {"attbigru2s", }:
        model = ModelAttRNN2s(args.seq_len, args.layer_rnn, args.class_num,
                              args.dropout_rate, args.hid_rnn,
                              args.n_vocab, args.n_embed,
                              is_stds=str2bool(args.is_stds),
                              model_type=args.model_type)
    elif args.model_type in {"transencoder", }:
        model = ModelTransEncoder(args.seq_len, args.layer_tfe, args.class_num,
                                  args.dropout_rate, args.d_model_tfe, args.nhead_tfe, args.nhid_tfe,
                                  args.n_vocab, args.n_embed,
                                  is_stds=str2bool(args.is_stds),
                                  model_type=args.model_type)
    elif args.model_type == "resnet18":
        model = ModelResNet18(args.class_num, args.dropout_rate, str2bool(args.is_stds))
    else:
        raise ValueError("model_type not right!")

    if use_cuda:
        model = model.cuda()

    if args.init_model is not None:
        print("loading pre-trained model: {}".format(args.init_model))
        para_dict = torch.load(args.init_model) if use_cuda else torch.load(args.init_model,
                                                                            map_location=torch.device('cpu'))
        model_dict = model.state_dict()
        model_dict.update(para_dict)
        model.load_state_dict(model_dict)

    # Loss and optimizer
    weight_rank = torch.from_numpy(np.array([1, args.pos_weight])).float()
    if use_cuda:
        weight_rank = weight_rank.cuda()
    criterion = nn.CrossEntropyLoss(weight=weight_rank)
    if args.optim_type == "Adam":
        optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    elif args.optim_type == "RMSprop":
        optimizer = torch.optim.RMSprop(model.parameters(), lr=args.lr)
    elif args.optim_type == "SGD":
        optimizer = torch.optim.SGD(model.parameters(), lr=args.lr, momentum=0.8)
    elif args.optim_type == "Ranger":
        # use Ranger optimizer
        # refer to https://github.com/lessw2020/Ranger-Deep-Learning-Optimizer
        # needs python>=3.6
        try:
            from utils.ranger2020 import Ranger
        except ImportError:
            raise ImportError("please check if ranger2020.py is in utils/ dir!")
        optimizer = Ranger(model.parameters(), lr=args.lr, betas=(0.95, 0.999), eps=1e-5)
    else:
        raise ValueError("optim_type is not right!")
    scheduler = StepLR(optimizer, step_size=args.lr_decay_step, gamma=args.lr_decay)

    # Train the model
    total_step = len(train_loader)
    print("total_step: {}".format(total_step))
    curr_best_accuracy = 0
    model.train()
    for epoch in range(args.max_epoch_num):
        curr_best_accuracy_epoch = 0
        tlosses = []
        start = time.time()
        for i, sfeatures in enumerate(train_loader):
            if args.model_type in {"bilstm", "bigru", "attbilstm", "attbigru",
                                   "transencoder", }:
                _, kmer, ipd_means, ipd_stds, pw_means, pw_stds, labels = sfeatures
                if use_cuda:
                    kmer = kmer.cuda()
                    ipd_means = ipd_means.cuda()
                    ipd_stds = ipd_stds.cuda()
                    pw_means = pw_means.cuda()
                    pw_stds = pw_stds.cuda()
                    labels = labels.cuda()
                # Forward pass
                outputs, logits = model(kmer, ipd_means, ipd_stds, pw_means, pw_stds)
                loss = criterion(outputs, labels)
                tlosses.append(loss.detach().item())
            elif args.model_type in {"resnet18", }:
                _, _, mats_ccs_mean, mats_ccs_std, labels = sfeatures
                if use_cuda:
                    mats_ccs_mean = mats_ccs_mean.cuda()
                    mats_ccs_std = mats_ccs_std.cuda()
                    labels = labels.cuda()
                # Forward pass
                outputs, logits = model(mats_ccs_mean, mats_ccs_std)
                loss = criterion(outputs, labels)
                tlosses.append(loss.detach().item())
            elif args.model_type in {"attbigru2s", }:
                _, kmer, ipd_means, ipd_stds, pw_means, pw_stds, kmer2, ipd_means2, ipd_stds2, pw_means2, pw_stds2, \
                    labels = sfeatures
                if use_cuda:
                    kmer = kmer.cuda()
                    ipd_means = ipd_means.cuda()
                    ipd_stds = ipd_stds.cuda()
                    pw_means = pw_means.cuda()
                    pw_stds = pw_stds.cuda()
                    kmer2 = kmer2.cuda()
                    ipd_means2 = ipd_means2.cuda()
                    ipd_stds2 = ipd_stds2.cuda()
                    pw_means2 = pw_means2.cuda()
                    pw_stds2 = pw_stds2.cuda()
                    labels = labels.cuda()
                # Forward pass
                outputs, logits = model(kmer, ipd_means, ipd_stds, pw_means, pw_stds,
                                        kmer2, ipd_means2, ipd_stds2, pw_means2, pw_stds2)
                loss = criterion(outputs, labels)
                tlosses.append(loss.detach().item())
            else:
                raise ValueError("model_type not right!")

            # Backward and optimize
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 0.5)
            optimizer.step()

            if (i + 1) % args.step_interval == 0 or i == total_step - 1:
                model.eval()
                with torch.no_grad():
                    vlosses, vlabels_total, vpredicted_total = [], [], []
                    for vi, vsfeatures in enumerate(valid_loader):
                        if args.model_type in {"bilstm", "bigru", "attbilstm", "attbigru",
                                               "transencoder", }:
                            _, vkmer, vipd_means, vipd_stds, vpw_means, vpw_stds, \
                                vlabels = vsfeatures
                            if use_cuda:
                                vkmer = vkmer.cuda()
                                vipd_means = vipd_means.cuda()
                                vipd_stds = vipd_stds.cuda()
                                vpw_means = vpw_means.cuda()
                                vpw_stds = vpw_stds.cuda()
                                vlabels = vlabels.cuda()
                            voutputs, vlogits = model(vkmer, vipd_means, vipd_stds, vpw_means, vpw_stds)
                            vloss = criterion(voutputs, vlabels)
                        elif args.model_type in {"attbigru2s", }:
                            _, vkmer, vipd_means, vipd_stds, vpw_means, vpw_stds, \
                                vkmer2, vipd_means2, vipd_stds2, vpw_means2, vpw_stds2, \
                                vlabels = vsfeatures
                            if use_cuda:
                                vkmer = vkmer.cuda()
                                vipd_means = vipd_means.cuda()
                                vipd_stds = vipd_stds.cuda()
                                vpw_means = vpw_means.cuda()
                                vpw_stds = vpw_stds.cuda()
                                vkmer2 = vkmer2.cuda()
                                vipd_means2 = vipd_means2.cuda()
                                vipd_stds2 = vipd_stds2.cuda()
                                vpw_means2 = vpw_means2.cuda()
                                vpw_stds2 = vpw_stds2.cuda()
                                vlabels = vlabels.cuda()
                            voutputs, vlogits = model(vkmer, vipd_means, vipd_stds, vpw_means, vpw_stds,
                                                      vkmer2, vipd_means2, vipd_stds2, vpw_means2, vpw_stds2)
                            vloss = criterion(voutputs, vlabels)
                        elif args.model_type in {"resnet18", }:
                            _, _, vmats_ccs_mean, vmats_ccs_std, vlabels = vsfeatures
                            if use_cuda:
                                vmats_ccs_mean = vmats_ccs_mean.cuda()
                                vmats_ccs_std = vmats_ccs_std.cuda()
                                vlabels = vlabels.cuda()
                            # Forward pass
                            voutputs, vlogits = model(vmats_ccs_mean, vmats_ccs_std)
                            vloss = criterion(voutputs, vlabels)
                        else:
                            raise ValueError("model_type not right!")

                        _, vpredicted = torch.max(vlogits.data, 1)

                        if use_cuda:
                            vlabels = vlabels.cpu()
                            vpredicted = vpredicted.cpu()
                        vlosses.append(vloss.item())
                        vlabels_total += vlabels.tolist()
                        vpredicted_total += vpredicted.tolist()

                    v_accuracy = metrics.accuracy_score(vlabels_total, vpredicted_total)
                    v_precision = metrics.precision_score(vlabels_total, vpredicted_total)
                    v_recall = metrics.recall_score(vlabels_total, vpredicted_total)
                    if v_accuracy > curr_best_accuracy_epoch:
                        curr_best_accuracy_epoch = v_accuracy
                        if curr_best_accuracy_epoch > curr_best_accuracy - 0.0002:
                            torch.save(model.state_dict(),
                                       model_dir + args.model_type + '.b{}_epoch{}.ckpt'.format(args.seq_len,
                                                                                                epoch + 1))

                    time_cost = time.time() - start
                    print('Epoch [{}/{}], Step [{}/{}], TrainLoss: {:.4f}; '
                          'ValidLoss: {:.4f}, '
                          'Accuracy: {:.4f}, Precision: {:.4f}, Recall: {:.4f}, '
                          'curr_epoch_best_accuracy: {:.4f}; Time: {:.2f}s'
                          .format(epoch + 1, args.max_epoch_num, i + 1, total_step, np.mean(tlosses),
                                  np.mean(vlosses), v_accuracy, v_precision, v_recall,
                                  curr_best_accuracy_epoch, time_cost))
                    tlosses = []
                    start = time.time()
                    sys.stdout.flush()
                model.train()
        scheduler.step()
        if curr_best_accuracy_epoch > curr_best_accuracy:
            curr_best_accuracy = curr_best_accuracy_epoch
        else:
            if epoch >= args.min_epoch_num - 1:
                print("best accuracy: {}, early stop!".format(curr_best_accuracy))
                break
    endtime = time.time()
    clear_linecache()
    print("[train]training cost {} seconds".format(endtime - total_start))


def main():
    parser = argparse.ArgumentParser("")
    st_input = parser.add_argument_group("INPUT")
    st_input.add_argument('--train_file', type=str, required=True)
    st_input.add_argument('--valid_file', type=str, required=True)

    st_output = parser.add_argument_group("OUTPUT")
    st_output.add_argument('--model_dir', type=str, required=True)

    st_train = parser.add_argument_group("TRAIN")
    # model param
    st_train.add_argument('--model_type', type=str, default="attbigru2s",
                          choices=["attbilstm", "attbigru", "bilstm", "bigru",
                                   "transencoder",
                                   "resnet18",
                                   "attbigru2s"],
                          required=False,
                          help="type of model to use, 'attbilstm', 'attbigru', "
                               "'bilstm', 'bigru', 'transencoder', 'resnet18', "
                               "'attbigru2s', "
                               "default: attbigru2s")
    st_train.add_argument('--seq_len', type=int, default=21, required=False,
                          help="len of kmer. default 21")
    st_train.add_argument('--is_stds', type=str, default="yes", required=False,
                          help="if using std features at ccs level, yes or no. default yes.")
    st_train.add_argument('--class_num', type=int, default=2, required=False)
    st_train.add_argument('--dropout_rate', type=float, default=0.5, required=False)

    # BiRNN/transformerencoder model param
    st_train.add_argument('--n_vocab', type=int, default=16, required=False,
                          help="base_seq vocab_size (15 base kinds from iupac)")
    st_train.add_argument('--n_embed', type=int, default=4, required=False,
                          help="base_seq embedding_size")

    # BiRNN model param
    st_train.add_argument('--layer_rnn', type=int, default=3,
                          required=False, help="BiRNN layer num, default 3")
    st_train.add_argument('--hid_rnn', type=int, default=256, required=False,
                          help="BiRNN hidden_size for combined feature")

    # transformerencoder model param
    st_train.add_argument('--layer_tfe', type=int, default=6,
                          required=False, help="transformer encoder layer num, default 6")
    st_train.add_argument('--d_model_tfe', type=int, default=256,
                          required=False, help="the number of expected features in the "
                                               "transformer encoder/decoder inputs")
    st_train.add_argument('--nhead_tfe', type=int, default=4,
                          required=False, help="the number of heads in the multiheadattention models")
    st_train.add_argument('--nhid_tfe', type=int, default=512,
                          required=False, help="the dimension of the feedforward network model")

    # model training
    st_train.add_argument('--optim_type', type=str, default="Adam", choices=["Adam", "RMSprop", "SGD",
                                                                             "Ranger"],
                          required=False, help="type of optimizer to use, 'Adam' or 'SGD' or 'RMSprop' or 'Ranger', "
                                               "default Adam")
    st_train.add_argument('--batch_size', type=int, default=512, required=False)
    st_train.add_argument('--lr', type=float, default=0.001, required=False)
    st_train.add_argument('--lr_decay', type=float, default=0.1, required=False)
    st_train.add_argument('--lr_decay_step', type=int, default=1, required=False)
    st_train.add_argument("--max_epoch_num", action="store", default=50, type=int,
                          required=False, help="max epoch num, default 50")
    st_train.add_argument("--min_epoch_num", action="store", default=10, type=int,
                          required=False, help="min epoch num, default 10")
    st_train.add_argument('--pos_weight', type=float, default=1.0, required=False)
    st_train.add_argument('--tseed', type=int, default=1234,
                          help='random seed for pytorch')
    st_train.add_argument('--step_interval', type=int, default=500, required=False)

    st_train.add_argument('--init_model', type=str, default=None, required=False,
                          help="file path of pre-trained model parameters to load before training")

    args = parser.parse_args()

    print("[main] start..")
    total_start = time.time()

    display_args(args)

    train(args)

    endtime = time.time()
    print("[main] costs {} seconds".format(endtime - total_start))


if __name__ == '__main__':
    main()
