import torch
from torch.utils.data import DataLoader
import torch.nn as nn
import torch.optim as optim
from torch.optim import lr_scheduler
import torchvision
from torchvision import transforms
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
from custom_dataset_research import ImagingDataset
from cifar_premodel import CIFAR, cifar10
from config import *

# TODO: load pretrained cifar-10 model, erase classifier, initialize classifier from scratch
# look at training code in altogether.py and follow the same patterns for ease of code
# build a config file for hyperparams

data_transforms = {
    'train': transforms.Compose([
        transforms.ToPILImage(),
        transforms.RandomHorizontalFlip(),
        transforms.ToTensor() #,
        # transforms.Normalize(mean=[0.485, 0.456, 0.406],
                             # std=[0.229, 0.224, 0.225])
    ]),
    'val': transforms.Compose([
        transforms.ToPILImage(),
        transforms.ToTensor() # ,
        # transforms.Normalize(mean=[0.485, 0.456, 0.406],
                             # std=[0.229, 0.224, 0.225])
    ]),
    'test': transforms.Compose([
        transforms.ToPILImage(),
        transforms.ToTensor() # ,
        # transforms.Normalize(mean=[0.485, 0.456, 0.406],
                             # std=[0.229, 0.224, 0.225])
    ]),
}
data_dir = './data/CRC/A02_180517_CellPatches/'
label_dir = './data/CRC/A03_180529_Labels/labels'

train_ids = [1919, 2019, 2075, 1843, 2036, 1915, 1932, 1928]
val_ids = [1794, 1939, 2038, 2382, 1798]
test_ids = [1792, 2107, 1910, 2023, 1879, 1957, 2090]

image_datasets = {x: ImagingDataset(data_dir, label_dir, train_ids, val_ids, test_ids, x, data_transforms[x]) for x in ['train', 'val']}
dataloaders = {x: DataLoader(image_datasets[x], batch_size = BATCH_SIZE, shuffle = True) for x in ['train', 'val']}
dataset_sizes = {x: len(image_datasets[x]) for x in ['train', 'val']}
print(dataset_sizes)
print(len(image_datasets['train'].images_dict))
print(torch.cuda.is_available())
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

loss_history = []
accuracy_history = []

def train_model(model, criterion, optimizer, scheduler, num_epochs = 50):
    since = time.time()

    best_model_wts = copy.deepcopy(model.state_dict())
    best_acc = 0.0
    for epoch in range(num_epochs):
        print('Epoch {}/{}'.format(epoch, num_epochs - 1))
        print('-' * 10)
        batch = 0
        for phase in ['train', 'val']:
            if phase == 'train':
                scheduler.step()
                model.train()
            else:
                model.eval()

            running_loss = 0.0
            running_corrects = 0
            num_batches = (1.0 * dataset_sizes[phase]) / BATCH_SIZE

            for inputs, labels in dataloaders[phase]:
                batch += 1
                # print(batch)
                if (batch % 100) == 0:
                    print('Batch ', batch, ' out of ', num_batches)
                inputs = inputs.to(device)
                labels = labels.to(device)

                optimizer.zero_grad()
                # t = time.time()
                with torch.set_grad_enabled(phase == 'train'):
                    outputs = model(inputs)
                    _, preds = torch.max(outputs, 1)
                    loss = criterion(outputs, labels)
                    if phase == 'train':
                        loss.backward()
                        optimizer.step()
                loss_history.append(loss.item())
                running_loss += loss.item() * inputs.size(0)
                running_corrects += torch.sum(preds == labels.data)
                # print(time.time() - t, 'forward-backward')
            epoch_loss = running_loss / dataset_sizes[phase]
            epoch_acc = running_corrects.double() / dataset_sizes[phase]

            print('{} Loss: {:.4f} Acc: {:.4f}'.format(
                phase, epoch_loss, epoch_acc))
            if phase == 'val':
                accuracy_history.append(epoch_acc)

            # deep copy the model
            if phase == 'val' and epoch_acc > best_acc:
                best_acc = epoch_acc
                best_model_wts = copy.deepcopy(model.state_dict())
        print()
    time_elapsed = time.time() - since
    print('Training complete in {:.0f}m {:.0f}s'.format(
        time_elapsed // 60, time_elapsed % 60))
    print('Best val Acc: {:4f}'.format(best_acc))

    # load best model weights
    model.load_state_dict(best_model_wts)
    return model

model_conv = cifar10(NUM_CHANNEL, None)
model_conv.classifier = None
model_conv.classifier = nn.Linear(8 * NUM_CHANNEL, 2)
nn.init.xavier_uniform_(model_conv.classifier.weight)
model_conv = model_conv.to(device)
criterion = nn.CrossEntropyLoss()
optimizer_conv = optim.SGD(model_conv.parameters(), lr = LEARNING_RATE, momentum = MOMENTUM)
exp_lr_scheduler = lr_scheduler.StepLR(optimizer_conv, step_size = STEP_SIZE, gamma = GAMMA)
model_conv = train_model(model_conv, criterion, optimizer_conv, exp_lr_scheduler, num_epochs = NUM_EPOCHS)

torch.save(model_conv, 'cnn_tumor_crc.pt')

np.savetxt('loss_history.txt', loss_history)
np.savetxt('acc.txt', accuracy_history)
torch.save(model_conv, 'conv')




