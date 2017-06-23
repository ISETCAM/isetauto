from PIL import Image
import numpy as np
import skimage.io as io
import tensorflow as tf
import xml.etree.ElementTree as ET
import os

FLAGS = tf.app.flags.FLAGS
tf.app.flags.DEFINE_string(
    'phase', None,
    'The name of phase - train/test')
tf.app.flags.DEFINE_string(
    'data_path', None,
    'Path of the data, e.g. /share/wandell/data/NN_Camera_Generalization/data/VOCdevkit/VOC2007/')
VOC_LABELS = {
    'none': (0, 'Background'),
    'aeroplane': (1, 'Vehicle'),
    'bicycle': (2, 'Vehicle'),
    'bird': (3, 'Animal'),
    'boat': (4, 'Vehicle'),
    'bottle': (5, 'Indoor'),
    'bus': (6, 'Vehicle'),
    'car': (7, 'Vehicle'),
    'cat': (8, 'Animal'),
    'chair': (9, 'Indoor'),
    'door': (9, 'Indoor'),
    'mat': (9, 'Indoor'),
    'cow': (10, 'Animal'),
    'diningtable': (11, 'Indoor'),
    'dog': (12, 'Animal'),
    'horse': (13, 'Animal'),
    'motorbike': (14, 'Vehicle'),
    'person': (15, 'Person'),
    'pottedplant': (16, 'Indoor'),
    'sheep': (17, 'Animal'),
    'sofa': (18, 'Indoor'),
    'train': (19, 'Vehicle'),
    'tvmonitor': (20, 'Indoor'),
}

def bytes_feature(value):
    if not isinstance(value, list):
        value = [value]
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=value))

def int64_feature(value):
    if not isinstance(value, list):
        value = [value]
    return tf.train.Feature(int64_list=tf.train.Int64List(value=value))

def float_feature(value):
    if not isinstance(value, list):
        value = [value]
    return tf.train.Feature(float_list=tf.train.FloatList(value=value))

def pascal_to_tfrecord(txtname,tfrecords_filename):
    # Input image/annotation path
    JPEGImages_file_path=FLAGS.data_path+'JPEGImages/'
    Annotations_file_path=FLAGS.data_path+'Annotations/'
    txtname=txtname
    tfrecords_filename=tfrecords_filename
    #Define a writer to write the tfrecord file
    writer = tf.python_io.TFRecordWriter(tfrecords_filename)

    with open(txtname, 'r') as f:
        image_index = [x.strip() for x in f.readlines()]
        
    for index in image_index:
        oldpath = os.path.join(FLAGS.data_path, 'JPEGImages', index + '.jpg')   
        #oldpath=os.path.join(JPEGImages_file_path,img_path)
        file_path = os.path.split(oldpath)
        lists = file_path[1].split('.')
        file_ext = lists[-1]
        file_pro = lists[-2]
        annotation_name=file_pro+'.xml'
        annotation_path=os.path.join(Annotations_file_path,annotation_name)
        img = np.array(Image.open(oldpath))
        print(annotation_path)
        image_data = tf.gfile.FastGFile(oldpath, 'r').read()
        bboxes = []
        xmin = []
        ymin = []
        xmax = []
        ymax = []
        difficult = []
        truncated = []
        labels_text = []
        labels = []
        filename =annotation_path
        tree = ET.parse(filename)
        root = tree.getroot()
        size = root.find('size')
        shape = [int(size.find('height').text),
                 int(size.find('width').text),
                 int(size.find('depth').text)]
        for obj in root.findall('object'):
            label = obj.find('name').text
            labels.append(int(VOC_LABELS[label][0]))
            labels_text.append(label.encode('ascii'))
            if obj.find('difficult'):
                difficult.append(int(obj.find('difficult').text))
            else:
                difficult.append(0)
            if obj.find('truncated'):
                truncated.append(int(obj.find('truncated').text))
            else:
                truncated.append(0)
            bbox = obj.find('bndbox')
            bboxes.append((float(bbox.find('ymin').text) / shape[0],
                           float(bbox.find('xmin').text) / shape[1],
                           float(bbox.find('ymax').text) / shape[0],
                           float(bbox.find('xmax').text) / shape[1]
                           ))
        for b in bboxes:
            assert len(b) == 4
            # pylint: disable=expression-not-assigned
            [l.append(point) for l, point in zip([ymin, xmin, ymax, xmax], b)]

        image_format = b'JPEG'
        example = tf.train.Example(features=tf.train.Features(feature={
                'image/height': int64_feature(shape[0]),
                'image/width': int64_feature(shape[1]),
                'image/channels': int64_feature(shape[2]),
                'image/shape': int64_feature(shape),
                'image/object/bbox/xmin': float_feature(xmin),
                'image/object/bbox/xmax': float_feature(xmax),
                'image/object/bbox/ymin': float_feature(ymin),
                'image/object/bbox/ymax': float_feature(ymax),
                'image/object/bbox/label': int64_feature(labels),
                'image/object/bbox/label_text': bytes_feature(labels_text),
                'image/object/bbox/difficult': int64_feature(difficult),
                'image/object/bbox/truncated': int64_feature(truncated),
                'image/format': bytes_feature(image_format),
                'image/encoded': bytes_feature(image_data)}))  
        writer.write(example.SerializeToString())
    writer.close()
    
def main(_):
    if (not FLAGS.phase) or (not FLAGS.data_path):
        raise ValueError('You must supply --phase and --data_path')

    # Generate tfrecord according to train/test phase
    if FLAGS.phase == 'train':
        txtname = os.path.join(FLAGS.data_path, 'ImageSets', 'Main',
                               'trainval.txt')
        tfrecords_filename = 'pascal_voc_'+'trainval'+'.tfrecords'
    else:
        txtname = os.path.join(FLAGS.data_path, 'ImageSets', 'Main',
                               'test.txt')
        tfrecords_filename = 'pascal_voc_'+'test'+'.tfrecords'
    pascal_to_tfrecord(txtname, tfrecords_filename)
    
if __name__ == '__main__':
    tf.app.run()

