import os

from aiogram import Bot, Dispatcher, types, executor

#модель
import torch
from PIL import Image
from RealESRGAN import RealESRGAN

def io(path_name):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    try:
        model = RealESRGAN(device, scale=4)
        model.load_weights('weights/RealESRGAN_x4.pth', download=True)

        image = Image.open(path_name).convert('RGB')

        sr_image = model.predict(image)

        #sr_image.save(f'results/sr_image_{random.random(1:10001)}.png')
        res_name = 'sr_image_'+path_name
        sr_image.save(res_name)

    except:
        return None
    return res_name



#бот

bot = Bot('6262803450:AAHJIi4vJyY8PcUubgK0gB1FVlcZFVF5DsU')
dp = Dispatcher(bot)

def firstkeyboard():
   button=[
      [types.KeyboardButton(text="/help"),
      types.KeyboardButton(text="/restart")]
   ]
   keyboard = types.ReplyKeyboardMarkup(
      keyboard=button,
      resize_keyboard=True)
   return keyboard



@dp.message_handler(commands=['start']) #Явно указываем в декораторе, на какую команду реагируем.
async def send_welcome(message: types.Message):
   await message.answer(f'Привет, {message.from_user.username}!\nЯ готов к работе)', reply_markup=firstkeyboard()) #Так как код работает асинхронно, то обязательно пишем await.
#-----------------------------------


@dp.message_handler(commands=['help'])
async def on_help_ask(message: types.Message):
    await bot.send_message(message.from_user.id, f'Я бот, исполюзующий нейросети для улучшения качества твоих '
                                                 f'фотографий. Загрузи сюда фотографию, качество которой хочешь '
                                                 f'улучшить, и я отправлю тебе картинку в более высоком разрешении.')


@dp.message_handler(commands=['restart'])
async def cancel_act(message: types.Message):
   if os.path.exists(f'{message.from_user.username}.png'):
      os.remove(f'{message.from_user.username}.png')
   await message.answer("Начнём заново.")
   await send_welcome(message)


@dp.message_handler(content_types=["photo", "document", "sticker"])
async def on_photo(message: types.Message):
    path_name = ''
    if message.photo:
        await bot.send_message(message.from_user.id, 'Принял')
        path_name = f'{message.from_user.username}.jpg'
        await message.photo[-1].download(path_name)
    elif message.document:
        await bot.send_message(message.from_user.id, 'Понял')
        path_name = f'{message.from_user.username}.jpg'
        await message.document.download(path_name)
    elif message.sticker:
        await bot.send_message(message.from_user.id, 'Осознал')
        path_name = f'{message.from_user.username}.jpg'
        await message.sticker.download(path_name)
    n_file = io(path_name)
    if n_file:
        await message.answer_document(open(n_file, "rb"))
        os.remove(f'C:\projects\RealHighResBot\{n_file}')
        os.remove(f'C:\projects\RealHighResBot\{path_name}')
    else:
        await message.reply("Сорри, слишком большая фотка. Я устал её обрабатывать, попробуй другую)")



@dp.message_handler()
async def echo(message: types.Message):
   await message.reply("""Извини, но я не понимаю этой команды! 
Воспользуйся командой /help для ознакомления с моим функционалом.""")

