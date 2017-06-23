import tensorflow as tf
import skimage.io as io
import matplotlib.pyplot as plt
slim = tf.contrib.slim

# Size of the input imnage
IMAGE_HEIGHT = 384
IMAGE_WIDTH = 384

features = {
    'image/encoded': tf.FixedLenFeature((), tf.string, default_value=''),
    'image/format': tf.FixedLenFeature((), tf.string, default_value='jpeg'),
    'image/height': tf.FixedLenFeature([1], tf.int64),
    'image/width': tf.FixedLenFeature([1], tf.int64),
    'image/channels': tf.FixedLenFeature([1], tf.int64),
    'image/shape': tf.FixedLenFeature([3], tf.int64),
    'image/object/bbox/xmin': tf.VarLenFeature(dtype=tf.float32),
    'image/object/bbox/ymin': tf.VarLenFeature(dtype=tf.float32),
    'image/object/bbox/xmax': tf.VarLenFeature(dtype=tf.float32),
    'image/object/bbox/ymax': tf.VarLenFeature(dtype=tf.float32),
    'image/object/bbox/label': tf.VarLenFeature(dtype=tf.int64),
    'image/object/bbox/difficult': tf.VarLenFeature(dtype=tf.int64),
    'image/object/bbox/truncated': tf.VarLenFeature(dtype=tf.int64),
}
items_to_handlers = {
    'image': slim.tfexample_decoder.Image('image/encoded', 'image/format'),
    'shape': slim.tfexample_decoder.Tensor('image/shape'),
    'object/bbox': slim.tfexample_decoder.BoundingBox(
            ['ymin', 'xmin', 'ymax', 'xmax'], 'image/object/bbox/'),
    'object/label': slim.tfexample_decoder.Tensor('image/object/bbox/label'),
    'object/difficult': slim.tfexample_decoder.Tensor('image/object/bbox/difficult'),
    'object/truncated': slim.tfexample_decoder.Tensor('image/object/bbox/truncated'),
}
decoder = slim.tfexample_decoder.TFExampleDecoder(features, items_to_handlers)

dataset = slim.dataset.Dataset(
        data_sources='pascal_voc_trainval.tfrecords',
        reader=tf.TFRecordReader,
        num_samples = 100,
        decoder=decoder,
        items_to_descriptions = {},
        num_classes=21)

provider = slim.dataset_data_provider.DatasetDataProvider(
                dataset,
                num_readers=100,
                shuffle=False)
[image, shape, glabels, gbboxes] = provider.get(['image', 'shape',
                                                             'object/label',
                                                             'object/bbox'])

# Random transformations can be put here: right before you crop images
# to predefined size. To get more information look at the stackoverflow
# question linked above.

resized_image = tf.image.resize_image_with_crop_or_pad(image=image,
                                       target_height=IMAGE_HEIGHT,
                                       target_width=IMAGE_WIDTH)


#shuffle_batch
images,shape= tf.train.shuffle_batch( [resized_image,shape],
                                             batch_size=10,
                                             num_threads=4,
                                             capacity=2000,
                                             min_after_dequeue=1000)
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

# The op for initializing the variables.
init_op = tf.group(tf.global_variables_initializer(),
                   tf.local_variables_initializer())

with tf.Session() as sess:
    
    sess.run(init_op)
    
    coord = tf.train.Coordinator()
    threads = tf.train.start_queue_runners(coord=coord)
    
    # Let's read off 3 batches just for example
    '''
    for i in range(3):
    
        img= sess.run([images,shape])
        img_batch = img[0]
        print(img_batch.shape)
        img_show=tf.reshape(img_batch[i, :, :, :], [IMAGE_HEIGHT, IMAGE_WIDTH,3])
        
        print('current batch')
        plt.imshow(sess.run(img_show), cmap='gray')
        plt.show()
    '''   
    for i in range(3):
    
        img,l= sess.run([images,shape])
        img_batch = img[0]
        print(img_batch.shape)
        print('current batch')
        print(l)
        plt.imshow(img_batch,cmap='gray')
        plt.show()
    coord.request_stop()
    coord.join(threads)
