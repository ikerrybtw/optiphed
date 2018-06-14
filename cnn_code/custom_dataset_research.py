from torch.utils.data import Dataset
from torchvision import transforms
import numpy as np
from skimage import io
import scipy.io
import os
import random
import time

class ImagingDataset(Dataset):
    def __init__(self, data_folder, label_folder, train_ids, val_ids, test_ids, mode='train', transform=None):
        self.train_ids = train_ids
        self.val_ids = val_ids
        self.test_ids = test_ids
        patient_folders = sorted(os.listdir(data_folder))
        self.train_dict = {}
        self.val_dict = {}
        self.test_dict = {}
        self.train_labels = {}
        self.val_labels = {}
        self.test_labels = {}
        self.train_keys = []
        self.val_keys = []
        self.test_keys = []
        # self.data_folder = data_folder
        # self.label_folder = label_folder
        for pat_id in patient_folders:
            print(pat_id)
            if int(pat_id) in self.train_ids:
                full_path = os.path.join(data_folder, pat_id)
                self.train_dict[int(pat_id)] = [os.path.join(full_path, image) for image in sorted(os.listdir(full_path))]
                full_label_path = os.path.join(label_folder, pat_id + '_labels.txt')
                f = open(full_label_path)
                pat_labels = f.readlines()
                pat_labels = [int(s.strip()) for s in pat_labels]
                self.train_labels[int(pat_id)] = pat_labels
                self.train_keys.append(int(pat_id))
            elif int(pat_id) in self.val_ids:
                full_path = os.path.join(data_folder, pat_id)
                self.val_dict[int(pat_id)] = [os.path.join(full_path, image) for image in sorted(os.listdir(full_path))]
                full_label_path = os.path.join(label_folder, pat_id + '_labels.txt')
                f = open(full_label_path)
                pat_labels = f.readlines()
                pat_labels = [int(s.strip()) for s in pat_labels]
                self.val_labels[int(pat_id)] = pat_labels
                self.val_keys.append(int(pat_id))
            elif int(pat_id) in self.test_ids:
                full_path = os.path.join(data_folder, pat_id)
                self.test_dict[int(pat_id)] = [os.path.join(full_path, image) for image in sorted(os.listdir(full_path))]
                full_label_path = os.path.join(label_folder, pat_id + '_labels.txt')
                f = open(full_label_path)
                pat_labels = f.readlines()
                pat_labels = [int(s.strip()) for s in pat_labels]
                self.test_labels[int(pat_id)] = pat_labels
                self.test_keys.append(int(pat_id))
            else:
                continue
        self.mode = mode
        if mode == 'train':
            self.images_dict = self.train_dict
            self.labels_dict = self.train_labels
            self.keys = self.train_keys
        elif mode == 'val':
            self.images_dict = self.val_dict
            self.labels_dict = self.val_labels
            self.keys = self.val_keys
        elif mode == 'test':
            self.images_dict = self.test_dict
            self.labels_dict = self.test_labels
            self.keys = self.test_keys
        else:
            raise NotImplementedError
        self.images_list = [f for folder in self.images_dict for f in self.images_dict[folder]]
        self.keys = sorted(self.keys)
        self.transform = transform
        self.lengths = []
        for key in self.keys:
            self.lengths.append(len(self.images_dict[key]))
        self.clengths = np.cumsum(self.lengths)

        # TODO: labels

    def __len__(self):
         return len(self.images_list)
         # return 50000 if self.mode == 'train' else 10000

    def __getitem__(self, index):
        tic = time.time()
        # for i in range(len(self.clengths)):
            # total = self.clengths[i]
            # if total > index:
                # break
        # key = self.keys[i]
        # i = random.randint(0, len(self.keys)-1)
        # key = self.keys[i]
        # ind = index - self.clengths[i-1] if i>0 else index
        # ind = random.randint(0, len(self.images_dict[key])-1)
        # print(key, self.lengths, ind, index)
        # image_list = self.images_dict[key]
        # image_name = image_list[random.randint(0, len(image_list)-1)]
        image_name = self.images_list[index] 
        str_lst = image_name.split('.')
        str_lst = str_lst[1].split('/')
        # print(str_lst)
        key = str_lst[-2]
        cell_id = str_lst[-1]
        # print(cell_id)
        # key = None
        image = io.imread(image_name)
        label = self.labels_dict[int(key)][int(cell_id)-1]
        if self.transform:
            image = self.transform(image)
        else:
            image = transforms.ToTensor(image)
        toc = time.time() - tic
        print(toc)
        return image, label





